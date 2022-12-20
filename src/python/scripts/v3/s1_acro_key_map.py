"""Generate acronym to key map

Usage:

python s1_acro_key_map.py \
    inp: data root
    inp: geojsons_path
    out: output path for acronym to key map

Example:

python src/python/scripts/v3/s1_acro_key_map.py \
    ~/Desktop/work/data/mouse_atlas \
    /v3/s0/geojson \
    /v3/s1/acro_key_map.json \

Created by Mukund on 2022-12-20
"""

import os
import sys
import csv
import json
import geopandas
from produtils import dprint

data_root = sys.argv[1]
geojsons_path = f'{data_root}{sys.argv[2]}'
acro_key_map_path = f'{data_root}{sys.argv[3]}'


files = os.listdir(geojsons_path)
files = [f for f in files if f.endswith('.geojson')]

acro_key_map = {}

for file in files:
    dprint(file)

    file_full_path = f'{geojsons_path}/{file}'

    with open(file_full_path) as f:
        data = json.load(f)

    for i in range(len(data["features"])):
        acro = data["features"][i]["properties"]["acronym"]
        key = data["features"][i]["id"]
        acro_key_map[acro] = key


# write acro key map as json
with open(acro_key_map_path, 'w') as f:
    json.dump(acro_key_map, f)

