"""
Generates (optionally) aggregated region labels for regionally aggregated and normalized gene
expression vis data generation in s9d_region_agg.py. Generates labels for all
CCF region as opposed to s9c_labelgen.py, which does a few selected regions.
Should be run before s9d_region_agg.py

Usage:

python s9d_labelgen.py
    inp: data root
    inp: path to folder with allen label ids
    out: path to folder to write out region labels

Usage example:

python src/python/scripts/analysis/s9d_labelgen.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels

Created by Mukund on 2022-07-07
"""

import sys
from produtils import dprint
from pathlib import Path
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.allen as allen
import csv

data_root = sys.argv[1]
ip_folder_labels = data_root+sys.argv[2]
op_folder = sys.argv[1]+sys.argv[3]


# get number of regions
for pid in range(1, 208, 2):

    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment


    # get region names/labels
    nis_id_str = str(apid).zfill(3)
    labels_csv_file = f'{ip_folder_labels}/allen_anno_data_{nis_id_str}.csv'
    dprint(labels_csv_file)
    region_ids = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_ids.append(int(row[9]))

    dprint(len(region_ids))
    aggr_labels_csv_file = f'{op_folder}/agg_labels_{nis_id_str}.csv'
    with open(aggr_labels_csv_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, row in enumerate(region_ids):
            line = [row]
            writer.writerow(line)
