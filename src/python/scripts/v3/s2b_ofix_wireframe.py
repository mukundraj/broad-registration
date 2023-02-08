"""Fixes pid 175-183 offset issue for new wireframes.

Usage:

python s2b_ofix_wireframe.py \
    io: data root
    ip: path to transformed wireframe folder
    op: path to ofixed wireframe folder


Example:

python src/python/scripts/v3/s2b_ofix_wireframe.py \
    ~/Desktop/work/data/mouse_atlas \
    /v3/s2/wireframes_improved_trans \
    /v3/s2/wireframes_ofixed \

Created by Mukund on 2023-02-07

"""
from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import csv
import numpy as np
import json
from src.python.utils.misc import get_ofixed_nis_idx
import os



data_root = sys.argv[1]
ip_folder_wireframes = f'{data_root}{sys.argv[2]}'
op_folder_wireframes = f'{data_root}{sys.argv[3]}'

# iterate over pids
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

    pid_str = str(pid).zfill(3)

    cnis_idx = str(get_ofixed_nis_idx(pid)).zfill(3)
    wireframe_file = f'{ip_folder_wireframes}/chuck_sp_wireframe_{cnis_idx}.png'

    op_wirekframe_file = f'{op_folder_wireframes}/chuck_sp_wireframe_{pid_str}.png'

    # copy wireframe to op_wirekframe_file
    cmd = f'cp {wireframe_file} {op_wirekframe_file}'
    os.system(cmd)



