"""Improves visibility of boundary of the wireframe

Usage:

python s2b_improve_wireframe.py
    inp: data root
    inp: path to folder with input wireframe files
    out: path to improved wireframes

Example:

python src/python/scripts/v3/s2b_improve_wireframe.py \
    ~/Desktop/work/data/mouse_atlas \
    /v3/s2/wireframes \
    /v3/s2/wireframes_improved \

Supplementary:

cd ~/Desktop/work/data/mouse_atlas/v3/s2/wireframes_improved
fd -x convert {} -transparent white ~/Desktop/work/data/mouse_atlas/v3/s2/wireframes_improved_trans/{}

Created by Mukund on 2023-01-03
"""


import sys
import os
from produtils import dprint
from PIL import Image, ImageFilter

wireframe_folder = f'{sys.argv[1]}{sys.argv[2]}'
out_folder = f'{sys.argv[1]}{sys.argv[3]}'

for file in os.listdir(wireframe_folder)[:]:
    if file.endswith(f'.png'):
        im_file = f'{wireframe_folder}/{file}'
        pid = int(file.split('.')[0].split('_')[-1])
        nis_id_str = str(pid).zfill(3)
        dprint(nis_id_str)

        # perform edge detection
        image = Image.open(im_file)

        # Find the edges by applying the filter ImageFilter.FIND_EDGES

        image = image.convert("L")

        imageWithEdges = image.filter(ImageFilter.FIND_EDGES)
        imageWireframe = imageWithEdges.point( lambda p: 0 if p > 0 else 255 )

        dprint(out_folder)
        wireframe_img_path = f'{out_folder}/chuck_sp_wireframe_{str(nis_id_str)}.png'
        # dprint(labelmap_img_path)
        # dprint(wireframe_img_path)
        imageWireframe.save(wireframe_img_path)



# cmd2 = f'gsutil cp {local_path} {bucket_path}'
# os.system(cmd2)
# dprint(cmd2)
