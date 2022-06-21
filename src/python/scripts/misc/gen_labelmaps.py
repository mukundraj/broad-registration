"""
Generates labelmap images in chuck space for Allen labels

Usage:

python src/python/scripts/misc/gen_labelmaps.py \
    inp: data root
    out: path to store labelmaps
    out: path to store wireframe

Usage example:

python src/python/scripts/misc/gen_labelmaps.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/allen_labels_imgs/labelmaps \
    /data_v3_nissl_post_qc/s9_analysis/allen_labels_imgs/wireframe \

python src/python/scripts/misc/gen_labelmaps.py \


Created by Mukund on 2022-06-15

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
from produtils import dprint
from PIL import Image, ImageFilter

data_root = sys.argv[1]
op_labelmaps = data_root+sys.argv[2]
op_wireframe = data_root+sys.argv[3]

id_to_color_map = allen.get_allen_regionid_to_color_map()

# dprint(id_to_color_map)

# create grid coords
axis2 = 4096
axis3 = 3606
resolution2 = 456
resolution3 = 320
x = np.linspace(0, axis2, resolution2*2)
y = np.linspace(0, axis3, resolution3*2)

xv, yv = np.meshgrid(x, y)

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

print("idx ", idx)

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
            # if (allen_ano > 0):
            #     return (0.1, 0.1, 0.9, 0.9)
            # else:
            #     return (0, 0, 0, 0.1)
            return id_to_color_map[allen_ano]
            # if (x > 630 and y>50):
            #     return (0.7, 0.7, 0.7, 0.9)
            # else:
            #     return (0, 0, 0, 0.1)

        return color_gen
    color_gen = color_gen_generator(allen_annos)

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

    labelmap_img_path = f'{op_labelmaps}/chuck_sp_labelmap_{str(nis_id_str)}.png'
    fig.savefig(labelmap_img_path, dpi=72, transparent=True)

    # Create an image object

    image = Image.open(labelmap_img_path)

    # Find the edges by applying the filter ImageFilter.FIND_EDGES

    image = image.convert("L")

    imageWithEdges = image.filter(ImageFilter.FIND_EDGES)
    imageWireframe = imageWithEdges.point( lambda p: 0 if p > 0 else 255 )

    wireframe_img_path = f'{op_wireframe}/chuck_sp_wireframe_{str(nis_id_str)}.png'
    # dprint(labelmap_img_path)
    # dprint(wireframe_img_path)
    imageWireframe.save(wireframe_img_path)
