import sys
import time
import threading

# === Spinner to show while app is loading ===
loading = True

def spinner():
    symbols = "|/-\\"
    idx = 0
    while loading:
        sys.stdout.write("\rLaunching GBIF Species Search... " + symbols[idx % len(symbols)])
        sys.stdout.flush()
        idx += 1
        time.sleep(0.1)
    sys.stdout.write("\rLaunching GBIF Species Search... done!        \n")

spinner_thread = threading.Thread(target=spinner)
spinner_thread.start()

# === Heavy imports happen here ===
from flask import Flask, request, render_template_string, send_file
from pygbif import occurrences, species as gbif_species
import math
import csv
import io
import os
from collections import OrderedDict
import concurrent.futures
import webbrowser
import threading as th

# === Stop spinner once imports are done ===
loading = False
spinner_thread.join()

# === The rest of your app code ===
app = Flask(__name__)

BASE_DIR = os.path.dirname(sys.executable if getattr(sys, 'frozen', False) else os.path.abspath(__file__))
CACHE_FILENAME = os.path.join(BASE_DIR, "species_cache.csv")

TAXA_GROUPS = {
    'Birds': 212,
    'Mammals': 359,
    'Amphibia': 131,
    'Lizards': 11592253,
    'Vascular plants': 7707728,
    'All plants': 6,
    'Fungi': 5,
    'Insects': 216,
    'Chordata': 44,
    'Animalia': 1,
    'Custom taxonKey': None
}

KINGDOM_SUBGROUPS = {
    1: [44,54,52,42,43,50,108,5967481,62,105],
    6: [35,13,7819616,36,106,9,7707728,220,194,7228684,245],
    5: [95,34,73,10791940,12238325,17,12252287,12237820,96,94,18],
}

HTML_TEMPLATE = """ 
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<title>GBIF Species Search</title>
<script>
function toggleCustomInput() {
    var taxonSelect = document.getElementById("taxon");
    var customInputDiv = document.getElementById("custom_taxon_div");
    if (taxonSelect.value === "Custom taxonKey") {
        customInputDiv.style.display = "block";
    } else {
        customInputDiv.style.display = "none";
    }
}
window.onload = toggleCustomInput;
</script>
<style>
  table {border-collapse: collapse; width: 100%;}
  th, td {border: 1px solid #ccc; padding: 6px; text-align: left;}
</style>
</head>
<body>
<h1>GBIF Species Search by Location</h1>
<form method="post">
    <label>Coordinates (lat, lon): 
        <input type="text" name="coordinates" required placeholder="-35.6264, 174.4823" value="{{coordinates or ''}}" />
    </label><br/>
    <label>Radius (km): <input type="number" name="radius" required min="1" max="50" value="{{radius or '5'}}" /></label><br/>
    <label>Taxonomic Group:
        <select id="taxon" name="taxon" onchange="toggleCustomInput()">
            {% for name in taxa_groups %}
                <option value="{{name}}" {% if name == taxon %}selected{% endif %}>{{name}}</option>
            {% endfor %}
        </select>
    </label><br/>
    <div id="custom_taxon_div" style="display:none;">
        <label>Enter custom taxonKey:
            <input type="text" name="custom_taxon" placeholder="e.g. 44" value="{{custom_taxon or ''}}" />
        </label><br/>
    </div>
    <input type="submit" value="Search" />
</form>

{% if error %}
<p style="color:red">{{error}}</p>
{% endif %}

{% if species_list %}
<h2>Species found ({{species_list|length}} unique):</h2>
<a href="/download">Download CSV</a>

<table>
    <thead>
        <tr>
            <th>Scientific Name</th>
            <th>Common Name</th>
            <th>Phylum</th>
            <th>Class</th>
            <th>Order</th>
            <th>Family</th>
        </tr>
    </thead>
    <tbody>
    {% for sp in species_list %}
        <tr>
            <td>{{ sp.scientific_name }}</td>
            <td>{{ sp.common_name }}</td>
            <td>{{ sp.phylum }}</td>
            <td>{{ sp.class_ }}</td>
            <td>{{ sp.order }}</td>
            <td>{{ sp.family }}</td>
        </tr>
    {% endfor %}
    </tbody>
</table>
{% endif %}
</body>
</html>
"""

def create_wkt_circle(lon, lat, radius_km, n_points=32):
    coords = []
    radius_deg = radius_km / 111  # approx degrees for latitude
    for i in range(n_points):
        angle = 2 * math.pi * i / n_points
        dx = radius_deg * math.cos(angle) / math.cos(math.radians(lat))
        dy = radius_deg * math.sin(angle)
        point_lon = lon + dx
        point_lat = lat + dy
        coords.append(f"{point_lon:.6f} {point_lat:.6f}")
    coords.append(coords[0])  # close polygon
    return "POLYGON((" + ",".join(coords) + "))"

species_cache = OrderedDict()
cache_dirty = False  # Track if cache has changed

def load_cache():
    global species_cache
    species_cache.clear()
    if not os.path.isfile(CACHE_FILENAME):
        return
    with open(CACHE_FILENAME, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row["speciesKey"]
            if key:
                species_cache[key] = {
                    "scientific_name": row["scientific_name"],
                    "common_name": row["common_name"],
                    "phylum": row["phylum"],
                    "class_": row["class"],
                    "order": row["order"],
                    "family": row["family"],
                }

def save_cache():
    with open(CACHE_FILENAME, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["speciesKey", "scientific_name", "common_name", "phylum", "class", "order", "family"])
        for key, data in species_cache.items():
            writer.writerow([
                key,
                data["scientific_name"],
                data["common_name"],
                data["phylum"],
                data["class_"],
                data["order"],
                data["family"],
            ])
    print(f"Cache saved to {CACHE_FILENAME}")

def append_cache_entry(species_key, data):
    global species_cache, cache_dirty
    if species_key not in species_cache:
        species_cache[species_key] = data
        cache_dirty = True

def evict_oldest_cache_entries(n):
    global species_cache
    if len(species_cache) <= n:
        return
    for _ in range(n):
        species_cache.popitem(last=False)
    save_cache()

def get_classification_and_common_name(species_key):
    try:
        data = gbif_species.name_usage(key=int(species_key))
        ranks = {"phylum": "", "class": "", "order": "", "family": ""}
        for rank in ranks.keys():
            if rank in data and data[rank]:
                ranks[rank] = data[rank]
        if not all(ranks.values()) and 'classification' in data:
            for rank_info in data['classification']:
                rank_name = rank_info.get('rank', '').lower()
                if rank_name in ranks and not ranks[rank_name]:
                    ranks[rank_name] = rank_info.get('name', '')
        common_name = data.get("vernacularName", "") or ""
        return ranks, common_name
    except Exception:
        return {"phylum": "", "class": "", "order": "", "family": ""}, ""

last_search_results = []

def query_species(lat, lon, radius_km, taxon_key):
    global last_search_results, cache_dirty

    if not taxon_key:
        return None, "No taxon key provided."

    wkt_polygon = create_wkt_circle(lon, lat, radius_km)
    keys_to_search = [taxon_key]
    if taxon_key in KINGDOM_SUBGROUPS:
        keys_to_search += KINGDOM_SUBGROUPS[taxon_key]

    species_dict = {}

    def fetch_for_key(t_key):
        try:
            response = occurrences.search(
                geometry=wkt_polygon,
                limit=300,
                offset=0,
                taxonKey=t_key
            )
            return response.get('results', [])
        except Exception:
            return []

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(fetch_for_key, k) for k in keys_to_search]
        for future in concurrent.futures.as_completed(futures):
            results = future.result()
            for rec in results:
                sci_name = rec.get('species')
                common_name = rec.get('vernacularName') or ""
                species_key = rec.get('speciesKey')
                if sci_name and sci_name not in species_dict:
                    species_dict[sci_name] = {
                        "common_name": common_name,
                        "species_key": species_key
                    }

    def enrich_species(item):
        sci, data = item
        species_key = data["species_key"]
        orig_common = data["common_name"]

        if species_key and species_key in species_cache:
            cached = species_cache[species_key]
            return {
                "scientific_name": cached["scientific_name"],
                "common_name": cached["common_name"],
                "phylum": cached["phylum"],
                "class_": cached["class_"],
                "order": cached["order"],
                "family": cached["family"],
            }

        if species_key:
            ranks, gbif_common = get_classification_and_common_name(species_key)
            common_name = gbif_common if gbif_common else orig_common
        else:
            ranks = {"phylum": "", "class": "", "order": "", "family": ""}
            common_name = orig_common

        data_to_cache = {
            "scientific_name": sci,
            "common_name": common_name,
            "phylum": ranks["phylum"],
            "class_": ranks["class"],
            "order": ranks["order"],
            "family": ranks["family"],
        }

        if species_key:
            append_cache_entry(species_key, data_to_cache)

        return data_to_cache

    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        enriched_species = list(executor.map(enrich_species, species_dict.items()))

    enriched_species.sort(key=lambda x: (
        x["phylum"] or "",
        x["class_"] or "",
        x["order"] or "",
        x["family"] or "",
        x["scientific_name"]
    ))

    # Save cache only if new entries added
    if cache_dirty:
        save_cache()
        cache_dirty = False

    # Trim cache only if exceeding 10000 entries
    MAX_CACHE_SIZE = 10000
    current_cache_size = len(species_cache)
    if current_cache_size > MAX_CACHE_SIZE:
        to_remove = current_cache_size - MAX_CACHE_SIZE
        evict_oldest_cache_entries(to_remove)

    last_search_results = enriched_species
    return enriched_species, None

@app.route("/", methods=["GET", "POST"])
def index():
    error = None
    species_list = []
    taxon = None
    custom_taxon = None
    radius = None
    coordinates = None

    if request.method == "POST":
        coordinates = request.form.get("coordinates", "").strip()
        radius = request.form.get("radius")
        taxon = request.form.get("taxon")
        custom_taxon = request.form.get("custom_taxon", "").strip()

        try:
            lat_str, lon_str = [s.strip() for s in coordinates.split(",")]
            latitude = float(lat_str)
            longitude = float(lon_str)
        except Exception:
            error = "Coordinates must be two numbers separated by a comma (e.g., -35.6264, 174.4823)."
            latitude = longitude = None

        try:
            radius = float(radius)
        except Exception:
            error = "Radius must be a number."

        taxon_key_to_use = None
        if taxon == "Custom taxonKey":
            if not custom_taxon.isdigit() or int(custom_taxon) <= 0:
                error = "Custom taxonKey must be a positive integer."
            else:
                taxon_key_to_use = int(custom_taxon)
        else:
            if taxon not in TAXA_GROUPS:
                error = "Invalid taxonomic group selected."
            else:
                taxon_key_to_use = TAXA_GROUPS[taxon]

        if not error:
            if not (-90 <= latitude <= 90) or not (-180 <= longitude <= 180):
                error = "Latitude must be between -90 and 90, longitude between -180 and 180."
            elif not (1 <= radius <= 50):
                error = "Radius must be between 1 and 50 km."
            else:
                species_list, error = query_species(latitude, longitude, radius, taxon_key_to_use)
                if species_list is None:
                    species_list = []
                    if error is None:
                        error = "No species found or API error."

    return render_template_string(
        HTML_TEMPLATE,
        species_list=species_list,
        error=error,
        taxon=taxon,
        custom_taxon=custom_taxon,
        radius=radius,
        coordinates=coordinates,
        taxa_groups=TAXA_GROUPS.keys()
    )

@app.route("/download")
def download_csv():
    global last_search_results
    if not last_search_results:
        return "No data to download. Please perform a search first.", 400

    si = io.StringIO()
    writer = csv.writer(si)
    writer.writerow(["Phylum", "Class", "Order", "Family", "Scientific Name", "Common Name"])
    for sp in last_search_results:
        writer.writerow([
            sp["phylum"],
            sp["class_"],
            sp["order"],
            sp["family"],
            sp["scientific_name"],
            sp["common_name"]
        ])

    output = io.BytesIO()
    output.write(si.getvalue().encode("utf-8"))
    output.seek(0)

    return send_file(
        output,
        mimetype="text/csv",
        download_name="results.csv",
        as_attachment=True,
    )

import webbrowser
import threading

def open_browser():
    webbrowser.open_new("http://127.0.0.1:5000")

if __name__ == "__main__":
    load_cache()
    threading.Timer(1, open_browser).start()
    # Turn off debug and reloader to avoid restart loops in PyInstaller
    app.run(debug=False, use_reloader=False)
