""" Script to render transformed bead positions with labeled beads in different
color.

Usage:

python input_csv_file output_png_file

python render_tfmed_poses.py \
    inp: nissl_id \
    inp: path to file mapping old and new nissl names (needed to determine labelmap nrrd filenames)
    inp: input_folder \
    out: output_folder 

Usage example:

python src/python/gists/render_tfmed_poses.py /Users/mraj/Desktop/sample-hz/output.csv /Users/mraj/Desktop/sample-hz/example.png

python src/python/gists/render_tfmed_poses.py \
    65 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s2_seg_ids/filenames_map.csv \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s4_bead_to_segid/bead_to_segid \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s4_bead_to_segid/bead_to_segid


Created updated version on 2022-04-22

"""

import cairo
import math
import numpy as np
import csv
import sys
from pathlib import Path
path_root = Path(__file__).parents[3]
sys.path.append(str(path_root))
import src.python.utils.io as io

def draw(surface, pts, radius, colors):
    for i, pt in enumerate(pts):
        x, y, radius = (pt[0], pt[1], radius)
        ctx = cairo.Context(surface)
        ctx.set_line_width(1)
        ctx.arc(x, y, radius, 0, 1.0 * math.pi)
        ctx.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
        ctx.fill_preserve()
        # ctx.set_source_rgb(1, 1, 1)
        ctx.set_source_rgb(colors[i][0], colors[i][1], colors[i][2])
        ctx.stroke()


nissl_id = int(sys.argv[1])
mapperfile_csv = sys.argv[2]
ip_folder = sys.argv[3]
op_folder = sys.argv[4]

img_dims = io.get_img_dimensions(mapperfile_csv)
nis_idx = str(nissl_id).zfill(3)

input_csv_file = f'{ip_folder}/bead_to_segid_{nis_idx}.csv'
op_file = f'{op_folder}/bead_to_segid_{nis_idx}.png'

print(input_csv_file)

# WIDTH, HEIGHT = 3253, 4643
WIDTH, HEIGHT = img_dims[nissl_id][0], img_dims[nissl_id][1]
# WIDTH, HEIGHT = 1728, 1728
# WIDTH, HEIGHT = 2048, 2048
# WIDTH, HEIGHT = 4096, 3606

print("WIDTH:", WIDTH, "HEIGHT:", HEIGHT)


surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
ctx = cairo.Context(surface)

# ctx.scale(WIDTH, HEIGHT)  # Normalizing the canvas

# pat = cairo.LinearGradient(0.0, 0.0, 0.0, 1.0)
# pat.add_color_stop_rgba(1, 0.7, 0, 0, 0.5)  # First stop, 50% opacity
# pat.add_color_stop_rgba(0, 0.9, 0.7, 0.2, 1)  # Last stop, 100% opacity

# ctx.rectangle(0, 0, 1, 1)  # Rectangle(x0, y0, x1, y1)
# ctx.set_source(pat)
# ctx.fill()

# ctx.translate(0.1, 0.1)  # Changing the current transformation matrix

# ctx.move_to(0, 0)
# # Arc(cx, cy, radius, start_angle, stop_angle)
# ctx.arc(0.2, 0.1, 0.1, -math.pi / 2, 0)
# ctx.line_to(0.5, 0.1)  # Line to (x,y)
# # Curve(x1, y1, x2, y2, x3, y3)
# ctx.curve_to(0.5, 0.2, 0.5, 0.4, 0.2, 0.8)
# ctx.close_path()

# ctx.set_source_rgb(0.3, 0.2, 0.5)  # Solid color
# ctx.set_line_width(0.02)
# ctx.stroke()

pts = np.random.rand(30, 2)*255
pts = list(pts)
colors = np.random.rand(30, 3)
radius = 5

# print(pts)

# draw(surface, pts, radius, colors)

# surface.write_to_png("example.png")  # Output to PNG

input_pts = []
labels = []
with open(input_csv_file, newline='\n') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # print(row['first_name'], row['last_name'])
        # row = list(map(float, row))
        point = [int(float(row[2])), int(float(row[3]))]
        input_pts.append(point)
        labels.append(int(row[0]))

# input_pts = np.array(input_pts)
# print(input_pts)
xmax, ymax = np.array(input_pts).max(axis=0)
print('xmax:', xmax,'ymax:', ymax)

colors = []
for label in labels:
    if (label==-1):
        colors.append([1,0,0])
    elif (label == 0):
        colors.append([0.5,0.5,0.5]) # unlabeled points
        # colors.append([0,1,0])
    else:
        colors.append([0, 1,0]) # labeled points

draw(surface, input_pts, radius, colors)

surface.write_to_png(op_file)  # Output to PNG
