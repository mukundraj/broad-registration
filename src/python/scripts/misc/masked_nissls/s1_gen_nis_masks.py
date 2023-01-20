"""Generate mask images for Nissl images to make region in Nissl images that lie outside the region sampled by slide-seq beads

Usage:

python s1_gen_nis_masks.py \
    io: data_root
    ip: coords in Chuck space
    op: mask images

Example:

python src/python/scripts/misc/masked_nissls/s1_gen_nis_masks.py \
    ~/Desktop/work/data/mouse_atlas \
    /v3/s2/coords \
    /misc/masked_nissls/s1/bead_masks \

Created by Mukund on 2023-01-17

"""

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import numpy as np
import sys
from produtils import dprint
import matplotlib.pyplot as plt
import cairo


io_data_root = sys.argv[1]
ip_coords_folder = f'{io_data_root}{sys.argv[2]}'
op_mask_folder = f'{io_data_root}{sys.argv[3]}'


def draw_polygon_with_cairo(pts, img_path):
    WIDTH, HEIGHT = 4096, 3605
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)
    ctx.save()
    ctx.set_source_rgba(0, 0, 0, 1)
    ctx.paint()
    ctx.restore()

    for i in range(len(pts)):
        if i == 0:
            ctx.move_to(pts[i][0], pts[i][1])
        else:
            ctx.line_to(pts[i][0], pts[i][1])
    ctx.close_path()
    ctx.set_source_rgba(1, 1, 1, 1)
    ctx.fill()
    surface.write_to_png(img_path)


# read in coords in Chuck space
start_pid = 1
end_pid = 207
pids = list(range(int(start_pid), int(end_pid)+1, 2))
if 5 in pids:
    pids.remove(5)
if 77 in pids:
    pids.remove(77)
if 167 in pids:
    pids.remove(167)
# iterate over pids
for pids_idx, pid in enumerate(pids):
    assert(pid!=5 and pid!=77 and pid!=167)

    coords_csv_name = f'{io_data_root}/v3/s2/coords/coords_{pid}.csv'
    dprint(f'coords_csv_name: {coords_csv_name}')

    # read csv file
    coords = []
    with (open(coords_csv_name, 'r')) as f:
        lines = f.readlines()
        for line in lines[1:]:
            line = line.split(':')[:-1]
            coords.append([int(line[0]), int(line[1])])
        coords = np.array(coords)
        dprint(coords)

    dprint(coords.shape)
    hull = ConvexHull(coords)

    # create convex hull image for mask
    # plt.plot(coords[hull.vertices,0], coords[hull.vertices,1], 'r--', lw=2)
    # plt.show()
    op_mask_img_path = f'{op_mask_folder}/mask_{pid}.png'
    draw_polygon_with_cairo(coords[hull.vertices], op_mask_img_path)



    # save convex hull image
