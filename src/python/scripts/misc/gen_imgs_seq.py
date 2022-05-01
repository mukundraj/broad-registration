""" Create a sequence of comparison images between allen CCF mouse brain and
beads transformed to CCF

Usage:

python gen_imgs_seq.py \
    out: path to grid pos folder \
    out: path to grid pos labels folder

Usage example:

python src/python/scripts/misc/gen_imgs_seq.py \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/chuck_sp_grid_pos \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/chuck_sp_grid_labels

Created by Mukund on 2022-04-27

Notes

matplotlib's imshow can take 

References

- https://stackoverflow.com/questions/66795706/how-can-i-iterate-over-a-meshgrid
- https://stackoverflow.com/questions/8218608/scipy-savefig-without-frames-axes-only-content
- https://stackoverflow.com/questions/52540037/create-image-using-matplotlib-imshow-meshgrid-and-custom-colors
- https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.axes.Axes.imshow.htmlkkkkkkkk

"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
from pathlib import Path
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.transforms as tfms
import src.python.utils.allen as allen

op_folder_grid = sys.argv[1]
op_folder_labels = sys.argv[2]

def density(x, y):
    return np.abs(x) + np.abs(y)

# create grid coords
axis2 = 4096
axis3 = 3606
resolution2 = 456
resolution3 = 320
x = np.linspace(0, axis2, resolution2*2)
y = np.linspace(0, axis3, resolution3*2)

xv, yv = np.meshgrid(x, y)

print(xv)

print(yv)


# def somefunc(x_value, y_value):
#     if x_value > 200:
#         return np.array(((1,0,0)))
#     elif y_value < 200:
#         return np.array(((0, 0, 1)))
#     else:
#         return np.array((0,1,0))

# res = np.zeros((xv.shape[0],xv.shape[1],3))
# for i in range(xv.shape[0]):
#     for j in range(xv.shape[1]):
#         res[i,j,:] = somefunc(xv[i,j],yv[i,j])

pts_grid = []
pos_map = {}
# write out positions to csv
idx = 0
for i in range(xv.shape[0]):
    for j in range(xv.shape[1]):
        # pts_grid.append([xv[i,j], yv[i,j]])
        pts_grid.append([xv[i,j], yv[i,j]])
        pos_map[(i,j)] = idx
        idx = idx + 1
        # print(i,j, xv[i,j], yv[i,j])

print("idx ",idx)


# op_file_grid = op_folder_grid+"/chuck_sp_grid_pts.csv"
# with open(op_file_grid, 'w', newline='\n') as csvfile:
#     writer = csv.writer(csvfile)
#     for idx, pt in enumerate(pts_grid):
#         row = [round(pt[0], 2), round(pt[1], 2)]
#         writer.writerow(row)


# transform to 3D space

print("before chuck to allend3d")
tfm_folder = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_fiducial_csvs" #fixme
tmp_storage = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/scratch" # fixme
qnii_json_file = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png_copy/QuickNii_4-17-22-mod.json" #fixme - hardcoding

nissl_ids = [i for i in np.arange(1,208,2)]
# nissl_ids = [i for i in np.arange(141,146,2)]
nissl_ids.remove(5)
nissl_ids.remove(77)
nissl_ids.remove(167)

# nissl_ids = [143]

for nissl_id in nissl_ids:

    nis_id_str = str(nissl_id).zfill(3)

    pts_allen = tfms.chuck_sp_to_allen3d(pts_grid, nissl_id, tfm_folder, qnii_json_file, tmp_storage)


    # sample back allen annotation id

    print("before sample allen annotation")
    nrrd_path = "/Users/mraj/Desktop/work/projects/active/broad-registration/annotation/ccf_2017/annotation_25.nrrd" # fixme - hardcoding
    allen_annos = allen.sample_allen_annotation(pts_allen, nrrd_path)
    print("after sample allen annotation")


    print(len(allen_annos), len(pts_grid))
    # create closure that uses annotation id to generate color

    def color_gen_generator(allen_annos):

        def color_gen(x, y):
            idx = pos_map[(x,y)]
            allen_ano = allen_annos[idx]
            if (allen_ano > 0):
                return (0.1, 0.1, 0.9, 0.9)
            else:
                return (0, 0, 0, 0.1)
            # if (x > 630 and y>50):
            #     return (0.7, 0.7, 0.7, 0.9)
            # else:
            #     return (0, 0, 0, 0.1)

        return color_gen
    color_gen = color_gen_generator(allen_annos)

    # print(xv)
    # print(np.shape(xv))
    # print("")
    # print(yv)

    # D = density(xv, yv)

    # print(D)

    print("xv shape", xv.shape)
    res = np.zeros((xv.shape[0],xv.shape[1],4))
    for i in range(xv.shape[0]):
        for j in range(xv.shape[1]):
            res[i,j,:] = color_gen(i,j)

    print("res shape", np.shape(res))

    # plt.pcolormesh(x, y, D)
    # plt.title('Density function = |x| + |y|')

    w=56.89 # from images of beads transformed to Chuck space
    h=50.08
    fig = plt.figure(frameon=False)
    fig.set_size_inches(w,h)
    # ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(res, aspect='auto')
    fig.savefig(f'{op_folder_labels}/chuck_sp_labelmap_{str(nis_id_str)}.png', dpi=72, transparent=True)


# transform grid coords to allen space

# sample region id and color for grid coords

# render rectangle block image for allen ccf

# plot bead positions as gray circles

# save





