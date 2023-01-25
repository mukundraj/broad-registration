""" Creates a packages zipped file containing nissls, and readme.txt to share
as on the docs.braincelldata.org/downloads page. The ppid_to_nissl.csv is just
for own reference.

Usage:

python s1a_nissls.py
    io: data root
    ip: id_to_tiff_mapper file
    ip: nissl src folder
    op: op_folder_name

Example:

python src/python/scripts/misc/for_downloads/s1a_nissls.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s2_seg_ids/id_to_tiff_mapper.csv \
    /data_v3_nissl_post_qc/s0_start_formatted_data/png_from_vsi \
    /misc/for_downloads/s1a_nissls \


Created by Mukund on 2023-01-25
"""

from produtils import dprint
import sys
import csv
import os
import shutil

from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import src.python.utils.filename_maps as fmaps



data_root = sys.argv[1]
id_to_tiff_mapper_file = f'{data_root}/{sys.argv[2]}'
nissl_src_folder = f'{data_root}/{sys.argv[3]}'
op_folder = f'{data_root}/{sys.argv[4]}'


pidToSrno, srnoToPid = fmaps.get_pidToSrno()

ppid_to_pid = {}

for pid, ppid in pidToSrno.items():
    dprint(pid, ppid)
    ppid_to_pid[ppid] = pid

pid_to_tiff_mapper = {}
# read csv id_to_tiff_mapper_file
with open(id_to_tiff_mapper_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        pid_to_tiff_mapper[row[0]] = row[1]


# create directory nissl in op_folder and overwrite if exists
nissl_folder = f'{op_folder}/nissls'
dprint(nissl_folder)
if os.path.exists(nissl_folder):
    shutil.rmtree(nissl_folder)
os.makedirs(nissl_folder)



for ppid in range(1, 102):
    pid = ppid_to_pid[str(ppid).zfill(3)]
    nis_file = f'{nissl_src_folder}/{pid_to_tiff_mapper[str(pid)]}'

    # replace extension of file name with .png
    nis_file_name = os.path.basename(nis_file)
    nis_file_name = nis_file_name.replace('.tiff', '.png')
    nis_file_name = f'nis_{nis_file_name}'

    nis_file_src = f'{nissl_src_folder}/{nis_file_name}'
    nis_file_dst = f'{nissl_folder}/{ppid}.png'

    # copy file from nis_file_src to nis_file_dst
    cmd = f'cp {nis_file_src} {nis_file_dst}'
    dprint(cmd)
    os.system(cmd)

# create a readme.txt file in op_folder
readme_file = f'{op_folder}/readme.txt'

with open(readme_file, 'w') as f:

    f.write("""The /nissls folder contains the nissl images seen in
    www.braincelldata.org in their original orientation. If you have any
    queries, please reach out at: braincelldata at gmail.com""")



# create a zip with containing nissl_folder and readme
cmd = f'cd {op_folder} && zip -r data_nissls.zip nissls readme.txt'
dprint(cmd)
os.system(cmd)


# write ppid to nissl file name mapping file
ppid_to_nissl_file_name = f'{op_folder}/ppid_to_nissl.csv'
with open(ppid_to_nissl_file_name, 'w') as f:
    for ppid in range(1, 102):
        pid = ppid_to_pid[str(ppid).zfill(3)]
        nis_file_name = os.path.basename(pid_to_tiff_mapper[str(pid)])
        nis_file_name = nis_file_name.replace('.tiff', '.png')
        nis_file_name = f'nis_{nis_file_name}'
        f.write(f'{ppid},{nis_file_name}\n')

