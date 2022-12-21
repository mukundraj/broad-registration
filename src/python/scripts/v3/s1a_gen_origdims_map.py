"""Generate csv file containing full resolution dimensions

Usage:

python s1_gen_origdims_map.py \
    inp: data root
    inp: id to tiff mapper file
    inp: vsi file path
    out: output path for csv file containing full resolution dimensions
    


Example:

python src/python/scripts/v3/s1_gen_origdims_map.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s2_seg_ids/id_to_tiff_mapper.csv \
    /data_v3_nissl_post_qc/s0_raw_data/vsi \
    /v3/s1/ \


Created by Mukund on 2022-12-21

"""

import os
import sys
import csv
from produtils import dprint


data_root = sys.argv[1]
id_to_tiff_file_path = f'{data_root}{sys.argv[2]}'
vsi_file_path = f'{data_root}{sys.argv[3]}'
orig_dims_map_path = f'{data_root}{sys.argv[4]}'
orig_dims_map_file = f'{orig_dims_map_path}/orig_dims_map.csv'

id_to_tiff = {}
with open(id_to_tiff_file_path) as f:
    reader = csv.reader(f)
    for row in reader:
        id_to_tiff[row[0]] = row[1]



dprint(f'orig_dims_map_path: {orig_dims_map_path}')

pids = list(id_to_tiff.keys())
dprint(pids)


orig_dims = []



for pid in pids:
    tiff_file = id_to_tiff[pid]
    vsi_filename = tiff_file.split(".")[0]+'.vsi'
    vsi_file = f'{vsi_file_path}/{vsi_filename}'
    dprint(f'vsi_file: {vsi_file}')
    cmd = f'/Users/mraj/Desktop/work/pkgs/bftools/showinf {vsi_file} > {orig_dims_map_path}temp/{pid}.txt'
    dprint(cmd)
    os.system(cmd)

    width = 0
    height = 0
    with open(f'{orig_dims_map_path}temp/{pid}.txt') as f:
        lines = f.readlines()
        for lno, line in enumerate(lines):

            if lno == 16:
                print(lno, line)
                width = line.split(" ")[-1].strip()
            elif lno == 17:
                print(lno, line)
                height = line.split(" ")[-1].strip()

        dprint(width, height)
        orig_dims.append([pid, vsi_filename, width, height])


with open(orig_dims_map_file, 'w') as f:
    writer = csv.writer(f)
    writer.writerows(orig_dims)









