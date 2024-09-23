"""
Plot required CCF region over points transformed to chuck space

python 04_plot_region.py \
    ip: region_name to plot
    ip: coords file with regions appended
    op: op img file path

python src/python/scripts/misc/for_jonahs/04_plot_region.py \
    'Dentate gyrus, granule cell layer' \
     "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2/tfmed_with_regions_048.csv" \
     "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2/tfmed_with_regions_048.tif" \

Created by Mukund on 2024-09-20
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

region_name = sys.argv[1]
coords_with_region_f = sys.argv[2]
op_img = sys.argv[3]

# read with | as delimiter 
coords_with_region = np.genfromtxt(coords_with_region_f, delimiter='|', dtype=None, encoding=None)

print('coords_with_region')
print(coords_with_region)


# get first two items in each tuple from coords_with_region if x[2].contains("Retrosp")==True

coords = np.array([(x[0], x[1]) for x in coords_with_region])
# get new extents
max_x = np.max(coords[:, 0])
max_y = np.max(coords[:, 1])
min_x = np.min(coords[:, 0])
min_y = np.min(coords[:, 1])
coords_reg = np.array([(x[0], x[1]) for x in coords_with_region if region_name in x[2] ])

print(coords)
print(np.shape(coords))

dpi=25.4

# plt.figure(figsize=(max_x-min_x, max_y-min_y))
plt.figure(figsize=(max_x/dpi, max_y/dpi))

# set fiture extents
plt.xlim(0, 4986)
plt.ylim(0, 3605)

plt.scatter(coords[:, 0], coords[:, 1], s=20, color='red')
plt.scatter(coords_reg[:, 0], coords_reg[:, 1], s=100, color='green')
plt.axis('off')

# make aspect ratio equal
plt.gca().set_aspect('equal', adjustable='box')

plt.savefig(op_img, bbox_inches='tight', pad_inches=0)


