"""
Prepares data full res nissl images and bead coordinates data for Josh Welch lab.


Example:

python src/python/scripts/misc/for_umich/prep_nissl_data.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s2_seg_ids/id_to_tiff_mapper.csv \
    /macosko_data/mraj/mouse_atlas_data/forJonah/tiff_from_vsi_maxres \
    /data_v3_nissl_post_qc/s4_bead_to_segid/bead_to_nis_coords \
    /misc/for_william \

Created by Mukund on 2023-04-11

"""

import sys
from produtils import dprint
import subprocess
import os
import shutil


data_root = sys.argv[1]
id_to_tiff_mapper_file = f'{data_root}/{sys.argv[2]}'
gbucket_tiff_folder = f'{sys.argv[3]}'
bead_to_nis_coords_folder = f'{data_root}/{sys.argv[4]}'
op_folder = f'{data_root}/{sys.argv[5]}'

# read id to tiff mapper csv
id_to_tiff_mapper = {}
with open(id_to_tiff_mapper_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line == '':
            continue
        seg_id, tiff_file = line.split(',')
        id_to_tiff_mapper[int(seg_id)] = tiff_file

# dprint(id_to_tiff_mapper)


# copy tiffs from gbucket to local
# iterate over pids
start_pid = 59 # 1
end_pid = 101 # 207
pids = list(range(int(start_pid), int(end_pid)+1, 2))

if 5 in pids:
    pids.remove(5)
if 77 in pids:
    pids.remove(77)
if 167 in pids:
    pids.remove(167)

# iterate over pids
for pids_idx, pid in enumerate(pids[:10]):

    pid_str = str(pid).zfill(3)
    tiff_file = id_to_tiff_mapper[pid]

    src = f'gs:/{gbucket_tiff_folder}/nis_{tiff_file}'
    dst = f'{op_folder}/pngs/{pid_str}.tif'
    download_cmd = f'gsutil cp {src} {dst}'
    convert_cmd = f'magick {dst} {op_folder}/pngs/{pid_str}.png'
    delete_cmd = f'rm {dst}'

    dprint(download_cmd)
    os.system(download_cmd)
    dprint(convert_cmd)
    os.system(convert_cmd)
    dprint(delete_cmd)
    os.system(delete_cmd)



dprint('done')

