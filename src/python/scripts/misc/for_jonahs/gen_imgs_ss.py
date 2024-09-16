"""

python src/python/scripts/misc/for_jonahs/gen_imgs_ss.py

"""

import numpy as np
import matplotlib.pyplot as plt
import subprocess

ip_file = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/ip_dlgn/dlgn_markers_bead_coords_df.csv"

# read csv
data = np.genfromtxt(ip_file, delimiter=',', dtype=None, encoding=None)

print(np.shape(data))

print(data)

coords = data[1:, 1:3]

# convert to a numeric
coords = coords.astype(np.float16)

colors = data[1:, 5]


print(coords)
print(colors)

opdir = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/op_dlgn"
op_img = f'{opdir}/ss.tif'
# plot scatterplot image with no padding



# get min/max of x and y coords
min_x = np.min(coords[:, 0])
max_x = np.max(coords[:, 0])

min_y = np.min(coords[:, 1])
max_y = np.max(coords[:, 1])


print(max_x, min_x, max_y, min_y)

# make coords positive
coords[:, 0] = coords[:, 0] - min_x
coords[:, 1] = coords[:, 1] - min_y

# get new extents
max_x = np.max(coords[:, 0])
max_y = np.max(coords[:, 1])
min_x = np.min(coords[:, 0])
min_y = np.min(coords[:, 1])

dpi=25.4
# plt.figure(figsize=(max_x-min_x, max_y-min_y))
plt.figure(figsize=(max_x/dpi, max_y/dpi))

# set fiture extents
plt.xlim(min_x, max_x)
plt.ylim(min_y, max_y)

plt.scatter(coords[:, 0], coords[:, 1], s=2500, color=colors)
plt.axis('off')

# make aspect ratio equal
plt.gca().set_aspect('equal', adjustable='box')

plt.savefig(op_img, bbox_inches='tight', pad_inches=0, dpi=dpi)
plt.close()

subprocess.run(["fd", "^.*tif$" , "-x", "magick", "{}", "-alpha", "off", "{}" ], cwd=opdir)
