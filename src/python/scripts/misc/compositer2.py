
"""
Script to plot aligned nissl over slice seq (to test alignment)

Usage:

python compositer2.py 
    inp: path to first folder \
    inp: path to second folder \
    out: path to output folder

Usage example:

python src/python/scripts/misc/compositer2.py \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/chuck_sp_grid_labels \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s3_registered_ss/chuck_space_recons_img \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/overlap_beads_on_ccf


"""
import sys
import subprocess
import numpy as np


nissl_path = sys.argv[1]
ss_path = sys.argv[2]
overlapped_path = sys.argv[3]

nissl_ids = [i for i in np.arange(1,208,2)]
nissl_ids.remove(5)
nissl_ids.remove(77)
nissl_ids.remove(167)

for id in nissl_ids:
    nis_idx = str(id).zfill(3)
    print(nis_idx)

    # if (id==175):
    #     corrected_nis_idx = str(173).zfill(3)
    # elif (id==177):
    #     corrected_nis_idx = str(175).zfill(3)
    # elif (id==179):
    #     corrected_nis_idx = str(177).zfill(3)
    # elif (id==181):
    #     corrected_nis_idx = str(179).zfill(3)
    # elif (id==183):
    #     corrected_nis_idx = str(181).zfill(3)
    # else:
    #     corrected_nis_idx = nis_idx

    nis_file = f'{nissl_path}/chuck_sp_labelmap_{nis_idx}.png'
    ss_file = f'{ss_path}/csp_recons_{nis_idx}.png'
    out_file = f'{overlapped_path}/overlapped_{nis_idx}.png'

    print(nis_file)
    subprocess.run([
        "composite",
        "-blend",
        "80",
        ss_file,
        nis_file,
        out_file
    ])

