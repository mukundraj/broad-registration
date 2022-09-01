"""
Script to generate a mapping from vsi filenames to final IDs in dataset

Usage:

python s2_fname_mapper_id_vsi.py
    inp: data root
    inp: path to filenames map file
    out: path to output file

Usage example:

python src/python/scripts/v2/s2_fname_mapper_id_vsi.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s2_seg_ids/filenames_map.csv \
    /data_v3_nissl_post_qc/s2_seg_ids/id_to_vsi_mapper.csv

Created by Mukund on 2022-08-30

"""
from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import csv

import src.python.utils.io as io
from produtils import dprint


data_root = sys.argv[1]
mapperfile_csv = data_root+sys.argv[2]
op_file = data_root+sys.argv[3]


mapper_to_new_filename, mapper_to_old_filename = io.get_filenames_map(mapperfile_csv)


keys = list(mapper_to_old_filename)

# op id_vsi_mapping.csv 

with open(op_file, 'w', newline='\n') as csvfile:
    writer = csv.writer(csvfile)
    for key in keys:
        filename = mapper_to_old_filename[key]
        filename = filename.replace("tif", "vsi").replace("nis_", "")
        row = [key, filename]
        writer.writerow(row)

