""" Copy the masked Nissl images to bucket for use in the portal.

Usage:

python s3_copy_maskd_nis.py \
    io: data_root
    ip: path to masked nissls
    op: path in bucket for GeneExp tab
    op: path in bucket for CellSpatial tab

Example:

python src/python/scripts/misc/masked_nissls/s3_copy_maskd_nis.py \
    ~/Desktop/work/data/mouse_atlas \
    /misc/masked_nissls/s2/masked_nissls \
    gs://bcdportaldata/genexp_data/gene_exprs_cshl \
    gs://bcdportaldata/cellspatial_data/cellscores \

Created by Mukund on 2023-01-18

"""


import sys
from produtils import dprint
import os

io_data_root = sys.argv[1]
ip_masked_nissls_folder = f'{io_data_root}{sys.argv[2]}'
op_geneexp_folder = sys.argv[3]
op_cellspatial_folder = sys.argv[4]


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

    local_path = f'{ip_masked_nissls_folder}/nis_{pid_str}.png'

    geneexp_path = f'{op_geneexp_folder}/puck{pid}/nis_{pid_str}.png'
    cmd = f'gsutil cp {local_path} {geneexp_path}'
    dprint(cmd)
    os.system(cmd)

    cellspatial_path = f'{op_cellspatial_folder}/puck{pid}/nis_{pid_str}.png'
    cmd = f'gsutil cp {local_path} {cellspatial_path}'
    dprint(cmd)
    os.system(cmd)

