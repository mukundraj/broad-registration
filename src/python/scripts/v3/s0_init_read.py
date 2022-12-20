"""Script to test plot registration data from Partha group.

Usage:

python init_read.py
    inp: data root
    inp: geojsons_path

Example:

python src/python/scripts/v3/init_read.py \
    ~/Desktop/work/data/mouse_atlas/v3 \
    /s0/geojson \

Created by Mukund Raj on 2022-12-19
"""

import geopandas
import matplotlib.pyplot as plt
import sys

data_root = sys.argv[1]
geojsons_path = f'{data_root}{sys.argv[2]}'



file = f'{geojsons_path}/atlas_to_MD743&742-N1-2019.02.13-21.15.52_MD743_1_0001.geojson'
# file = '/Users/mraj/Desktop/debug/sl1.geojson'
print(file)

data = geopandas.read_file(file)
print(data.head())
data.head()
data.plot(aspect=1) # https://stackoverflow.com/questions/67693007/valueerror-box-aspect-and-fig-aspect-must-be-positive
plt.show()

