
"""
Prepares bead coordinates data for Josh Welch lab. Also transforms bead coords to orig nissl space.


Example:

python src/python/scripts/misc/for_umich/prep_bead_coords.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s4_bead_to_segid/bead_to_nis_coords \
    /misc/for_william/pngs \
    /misc/for_william/bead_to_nis_coords_fullres \

Created by Mukund on 2023-04-18

"""

import sys
from produtils import dprint
import subprocess
import os
import shutil
import numpy as np
import PIL.Image
from PIL import Image
import matplotlib.pyplot as plt

data_root = sys.argv[1]
inp_bead_coords_folder = f'{data_root}/{sys.argv[2]}'
pngs_folder = f'{data_root}/{sys.argv[3]}'
op_folder = f'{data_root}/{sys.argv[4]}'

# copy bead to nis coords folder
# cmd = f'cp -r {bead_to_nis_coords_folder} {op_folder}'
# dprint(cmd)
# os.system(cmd)

start_pid = 1
end_pid = 207 # 207
pids = list(range(int(start_pid), int(end_pid)+1, 2))

if 5 in pids:
    pids.remove(5)
if 77 in pids:
    pids.remove(77)
if 167 in pids:
    pids.remove(167)

# iterate over pids
for pids_idx, pid in enumerate(pids):

    pid_str = str(pid).zfill(3)
    # read csv
    ip_csv = f'{inp_bead_coords_folder}/bead_to_nis_coords_{pid_str}.csv'

    # read csv using numpy skip first column but keep header
    data = np.genfromtxt(ip_csv, delimiter=',', skip_header=0, usecols=range(1, 5), dtype=np.int32)

    PIL.Image.MAX_IMAGE_PIXELS = 2*492324261 # a very large number

    # read full res image sizes
    img_file = f'{pngs_folder}/{pid_str}.png'
    img = Image.open(img_file)
    width, height = img.size


    coords_full = np.zeros((data.shape[0], 2), dtype=np.int32)
    # convert bead coords to orig nissl space
    coords_full[:,0] = width*(data[:,0]/data[:,2])
    coords_full[:,1] = height*(data[:,1]/data[:,3])

    # typecase coords_full to int
    coords_full = coords_full.astype(np.int32)

    # write csv
    op_csv = f'{op_folder}/bead_to_nis_coords_{pid_str}.csv'
    np.savetxt(op_csv, coords_full, delimiter=',', fmt='%d', header='x,y', comments='')

    # pick 20 percent of indices    
    idxs = np.random.choice(coords_full.shape[0], int(0.2*coords_full.shape[0]), replace=False)

    plt.figure()
    # plot 20pct beads for sanity check
    plt.scatter(coords_full[idxs,0], height-coords_full[idxs,1], s=0.5, alpha=0.5)
    plt.savefig(f'{op_folder}/bead_to_nis_coords_{pid_str}.png')

