"""
Script to render beads in Chuck space

Created by Mukund on 2022-03-22

Usage:

python render_beads_chuck_space.py \
    nissl_id or negative value for processing list of ids
    width height
    path to folder with positions
    path to output folder to store imgs

Usage example:

python src/python/gists/render_beads_chuck_space.py \
    143 \
    4096 3606 \
    /Users/mraj/Desktop/forgif/chuck_space_img_coords \
    /Users/mraj/Desktop/forgif/chuck_space_recons_img \

"""
from pathlib import Path
import sys
path_root = Path(__file__).parents[3]
sys.path.append(str(path_root))

import cairo
import src.python.utils.render2d as render2d
import sys
import csv
import numpy as np

nissl_id = int(sys.argv[1])
WIDTH = int(sys.argv[2])
HEIGHT = int(sys.argv[3])
ip_folder = sys.argv[4]
op_folder = sys.argv[5]


# WIDTH, HEIGHT = 3253, 4643


if (nissl_id > -1):
    nissl_ids = [nissl_id]
else:
    # nissl_ids = [123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143]
    nissl_ids = [133, 135, 137, 139, 143]
    nissl_ids = [1, 3, 133, 135, 137, 139, 143]
    nissl_ids = [i for i in np.arange(1,228,2)]
    nissl_ids.remove(5)
    nissl_ids.remove(77)
    nissl_ids.remove(167)

for nissl_id in nissl_ids:
    nis_id_str = str(nissl_id).zfill(3)
    ip_file = ip_folder+"/chuck_sp_img_coords_"+nis_id_str+".csv"
    print(f"Reading file {ip_file}")

    # read input coords
    input_pts = []
    with open(ip_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(int, row))
            point = [row[0], row[1]]
            input_pts.append(point)

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)
    colors =[[0.5, 0.5, 0.5]]
    radius = 2
    render2d.draw_disk(surface, input_pts, radius, colors)

    op_file = op_folder+"/csp_recons_"+nis_id_str+".png"
    surface.write_to_png(op_file)  # Output to PNG
