"""
Script for creating symlinks for viewer2 data

Usage:

python symlink_creator.py \
    symlink start path \
    symlink end path

Usage example:

python src/python/scripts/misc/symlink_creator.py \
    /Users/mraj/Desktop/work/projects/active/broad-registration/src/viewer2/data \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_allen_coords


Created by Mukund on 2022-04-10

"""

import sys
import numpy as np
import subprocess



nissl_ids = [i for i in np.arange(1,228,2)]

source_path = sys.argv[1]
target_path = sys.argv[2]

for nis_id in nissl_ids:
    nis_id_str = str(nis_id).zfill(3)
    source = f'./allen_img_coords_{nis_id_str}.csv'
    target = f'{target_path}/allen_img_coords_{nis_id_str}.csv'
    print(source)
    subprocess.run(["ln", "-s", target, source], cwd=source_path)


# subprocess.run(["fd", "^.*tif$" , "-x", "convert", "{}", "-alpha", "off", "{}" ], cwd=op_path)


