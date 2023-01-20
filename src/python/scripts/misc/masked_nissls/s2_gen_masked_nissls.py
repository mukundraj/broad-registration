"""

python src/python/scripts/misc/masked_nissls/s2_gen_masked_nissls.py \
    ~/Desktop/work/data/mouse_atlas \
    /misc/masked_nissls/s1/bead_masks \
    /misc/masked_nissls/s1/ccf_masks \
    /data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png_ofix \
    /misc/masked_nissls/s2/temp \
    /misc/masked_nissls/s2/masked_nissls \

Created by Mukund on 2023-01-18

"""

import sys
import os
from produtils import dprint

io_data_root = sys.argv[1]
ip_bead_masks = f'{io_data_root}{sys.argv[2]}'
ip_ccf_masks = f'{io_data_root}{sys.argv[3]}'
ip_ofix_nissls = f'{io_data_root}{sys.argv[4]}'
io_temp = f'{io_data_root}{sys.argv[5]}'
op_masked_nissls = f'{io_data_root}{sys.argv[6]}'


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
    assert(pid!=5 and pid!=77 and pid!=167)

    dprint(pid)
    pid_str = str(pid).zfill(3)

    # generate intersection mask of bead region and Allen CCF region
    bead_mask = f'{ip_bead_masks}/mask_{pid}.png'
    ccf_mask = f'{ip_ccf_masks}/chuck_sp_labelmap_{pid_str}.png'
    tmp_intersection_file = f'{io_temp}/tmp_bead_mask.png'
    cmd = f'composite -compose Multiply {bead_mask} {ccf_mask} {tmp_intersection_file}'
    os.system(cmd)

    # generate a hard masked nissl
    ofix_nissl = f'{ip_ofix_nissls}/nis_{pid_str}.png'
    hard_masked_nissl = f'{io_temp}/tmp_hard_masked_nissl.png'
    cmd = f'convert {ofix_nissl} {tmp_intersection_file} -alpha off -compose copy-opacity -composite {hard_masked_nissl}'
    os.system(cmd)

    # generate a soft masked nissl by blending ofixed nissl and hard masked nissl
    masked_nissl = f'{op_masked_nissls}/nis_{pid_str}.png'
    cmd = f'magick {hard_masked_nissl} {ofix_nissl} -define compose:args=34,76 -compose blend -composite {masked_nissl}'
    os.system(cmd)
