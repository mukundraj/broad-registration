"""
Script to update names of beadcoords files to match that on bcd portal.

python src/python/scripts/misc/for_jonahs/01b_get_all_beads_chksp.py

Created by Mukund on 2024-09-23
"""

import sys
from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import os
import subprocess
from pathlib import Path
import src.python.utils.filename_maps as fm



ipdir = Path("~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds").expanduser()

opdir = Path("/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/chuck_img_coords_allbds").expanduser()


pidToSrno, srnoToPid = fm.get_pidToSrno()


# get all pids
pids = pidToSrno.keys()

print (len(pids))

# iterater over all pids
for pid in pids:
    pid_padded = str(pid).zfill(3)
    print (pid, pid_padded)

    from_file = f'{ipdir}/chuck_sp_img_coords_{pid_padded}.csv'
    to_file = f'{opdir}/chuck_sp_img_coords_{pidToSrno[pid]}.csv'

    # copy from to to
    print (f'cp {from_file} {to_file}')

    os.system(f'cp {from_file} {to_file}')

