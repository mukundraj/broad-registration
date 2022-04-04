"""
Script to plot aligned nissl over slice seq (to test alignment)

Usage example:

python src/python/scripts/misc/compositer.py /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png /Users/mraj/Desktop/forgif/chuck_space_recons_img /Users/mraj/Desktop/forgif/overlapped

Created by Mukund on 2022-04-04

"""
import sys
import subprocess


nissl_path = sys.argv[1]
ss_path = sys.argv[2]
overlapped_path = sys.argv[3]

ids = list(range(1, 228, 2))

for id in ids:
    nis_idx = str(id).zfill(3)
    print(nis_idx)

    nis_file = f'{nissl_path}/nis_{nis_idx}.png'
    ss_file = f'{ss_path}/csp_recons_{nis_idx}.png'
    out_file = f'{overlapped_path}/overlapped_{nis_idx}.png'

    print(nis_file)
    subprocess.run([
        "composite",
        ss_file,
        nis_file,
        out_file
    ])

