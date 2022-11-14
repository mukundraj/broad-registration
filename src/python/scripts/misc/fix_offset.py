"""

Creates set of renamed images (nissl and wireframe) that takes into account the mapping offset issue in puckids 175-183.

Usage:

python fix_offset.py
    i/o: data_root
    inp: path to nissl images folder
    inp: path to atlas wireframe images folder
    out: path to renamed nissl images folder
    out: path to renamed atlas wireframe images folder

Usage example:

python src/python/scripts/misc/fix_offset.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png \
    /data_v3_nissl_post_qc/s7_annotations/allen_labels_imgs/wireframe_trans_bg \
    /data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png_ofix \
    /data_v3_nissl_post_qc/s7_annotations/allen_labels_imgs/wireframe_trans_bg_ofix \


Created by Mukund on 2022-11-11

"""

from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import src.python.utils.misc as misc
import subprocess

from produtils import dprint


data_root = sys.argv[1]
ip_nis = sys.argv[1] + sys.argv[2]
ip_wireframe = sys.argv[1] + sys.argv[3]
op_nis = sys.argv[1] + sys.argv[4]
op_wireframe = sys.argv[1] + sys.argv[5]



ids = list(range(1, 208, 2))

for id in ids:
    cnis_idx = str(misc.get_ofixed_nis_idx(id)).zfill(3)
    dprint(id, cnis_idx)
    nis_cmd = ['cp', ip_nis + '/nis_' + str(id).zfill(3) + '.png', op_nis + '/nis_' + cnis_idx + '.png']
    wireframe_cmd = ['cp', ip_wireframe + '/chuck_sp_wireframe_' + str(id).zfill(3) + '.png', op_wireframe + '/chuck_sp_wireframe_' + cnis_idx + '.png']

    subprocess.run(nis_cmd)
    subprocess.run(wireframe_cmd)

