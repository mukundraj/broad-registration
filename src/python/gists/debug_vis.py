"""
Script to visualize coordinates from CSV

Created by Mukund on 2022-03-23

Usage:

python debug_vis.py \
    nissl_id or negative value for processing list of ids
    width height
    path to folder with positions
    path to output folder to store imgs

Usage example:

python src/python/gists/debug_vis.py \
    135 \
    4000 4000 \
    /Users/mraj/Desktop/forgif/chuck_space_img_coords \
    /Users/mraj/Desktop/forgif/chuck_space_recons_img \

python src/python/gists/debug_vis.py \
    45 \
    2048 2048 \
    /Users/mraj/Desktop/forgif/fordebug \
    /Users/mraj/Desktop/forgif/fordebug \

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

for nissl_id in nissl_ids:
    nis_id_str = str(nissl_id).zfill(3)
    ip_file = ip_folder+"/chuck_sp_img_coords_"+nis_id_str+".csv"
    print(f"Reading file {ip_file}")

    # read input coords
    input_pts = []
    with open(ip_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            row = list(map(float, row))
            # row = list(map(lambda val: int((2+val)*1000), row)) # pts_tfmed_ss
            # row = list(map(lambda val: int(-2048*(val)), row))
            row = list(map(lambda val: int((val)), row))
            point = [row[0], row[1]]
            print(point)
            input_pts.append(point)

    print(np.amax(np.array(input_pts)), np.amin(np.array(input_pts)))
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)
    colors =[[0.5, 0.5, 0.5]]
    radius = 2
    render2d.draw_disk(surface, input_pts, radius, colors)

    op_file = op_folder+"/csp_recons_"+nis_id_str+".png"
    surface.write_to_png(op_file)  # Output to PNG
