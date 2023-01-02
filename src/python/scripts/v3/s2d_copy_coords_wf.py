"""Copy wireframe from local folder to puck specific folders in bucket (after
backing up old version). Copy to both gene and celltype related folders.

Usage:

python s2d_copy_coords_wf.py
    inp: data_root
    inp: path to new wireframes
    inp: path to new coords
    out: 
Example:

python src/python/scripts/v3/s2d_copy_coords_wf.py \
    ~/Desktop/work/data/mouse_atlas \
    /v3/s2/wireframes_trans \
    /v3/s2/coords   \
    gs://bcdportaldata/genexp_data/gene_exprs \
    /v3/s2/backup \

Created by Mukund on 2022-12-21
"""

from produtils import dprint
import os
import sys


data_root = sys.argv[1]
ip_folder_wireframes = f'{data_root}{sys.argv[2]}'
ip_folder_coords = f'{data_root}{sys.argv[3]}'
op_genexp_folder = sys.argv[4]
backup_dir = f'{data_root}{sys.argv[5]}'



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

    wireframe_file = f'chuck_sp_wireframe_{pid_str}.png'
    bucket_path = f'{op_genexp_folder}/puck{pid}/{wireframe_file}'
    bucket_path_coords = f'{op_genexp_folder}/puck{pid}/coords.csv'

    # # backup wireframe
    # backup_localpath = f'{backup_dir}/wireframes/{wireframe_file}'
    # cmd1 = f'gsutil -m cp {bucket_path} {backup_localpath}'
    # os.system(cmd1)

    # # backup coords
    # coords_file = f'coords_{pid}.csv'
    # backup_localpath = f'{backup_dir}/coords/{coords_file}'
    # cmd3 = f'gsutil -m rsync {bucket_path_coords} {backup_localpath}'
    # dprint(cmd3)
    # os.system(cmd3)


    # copy wireframe to gene related folder
    local_path = f'{ip_folder_wireframes}/{wireframe_file}'
    cmd2 = f'gsutil cp {local_path} {bucket_path}'
    os.system(cmd2)
    dprint(cmd2)

    # copy coords file to gene related folder
    local_path_coords = f'{ip_folder_coords}/coords_{pid}.csv'
    cmd4 = f'gsutil cp {local_path_coords} {bucket_path_coords}'
    dprint(cmd4)
    os.system(cmd4)
