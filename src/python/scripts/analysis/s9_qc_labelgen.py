"""
Generates aggregated labels for qc heatmap visualization. Starting by 
aggregating cortex layer labels. Run before s9c_qc_main.py

Usage:

python s9_qc_labelgen.py \
    inp: data root
    inp: path to folder with allen label ids \
    out: path to folder to write out aggregated labels

Usage example:

python src/python/scripts/analysis/s9_qc_labelgen.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s9_analysis/aggregated_labels

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

regions_ids = allen.get_cortex_layer_and_hippo_ids_lists()
dprint(len(regions_ids))
dprint(regions_ids[0])


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
    # [0 if x in [1,2] else 1 for x in [1,2,3]]
    region_ids = [-1 if x in regions_ids[0] else x for x in region_ids]
    region_ids = [-2 if x in regions_ids[1] else x for x in region_ids]
    region_ids = [-3 if x in regions_ids[2] else x for x in region_ids]
    region_ids = [-4 if x in regions_ids[3] else x for x in region_ids]
    region_ids = [-5 if x in regions_ids[4] else x for x in region_ids]
    region_ids = [-6 if x in regions_ids[5] else x for x in region_ids]

    region_ids = [-x if (int(x) < 0) else 99999 for x in region_ids]

    # non = [x for x in region_ids if x != 99999]
    # dprint(region_ids)

    aggr_labels_csv_file = f'{op_folder}/agg_labels_{nis_id_str}.csv'
    with open(aggr_labels_csv_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, row in enumerate(region_ids):
            line = [row]
            writer.writerow(line)


