
"""
Script for processing all in-tissue beads across all available pucks and
preparing a dict of region names mapping to position in list.

Usage:

python s9a_region_names.py  \
    inp: folder path to bead_ccf_labels_allbds with bead metadata in csv format
    inp: folder path to integrated_mats folder with processed gene counts in annodata h5ad format
    out: file path for writing labels dict containing gene names and index

Usage example:

python src/python/scripts/analysis/s9a_region_names.py \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/ccf_regions.json

Created by Mukund on 2022-05-19

References - https://stackoverflow.com/questions/9001509/how-can-i-sort-a-dictionary-by-key

"""

import sys
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix, hstack
import zarr
import json
from produtils import dprint
import os
import shutil
import csv
from collections import OrderedDict

ip_folder_labels = sys.argv[1]
ip_folder_counts = sys.argv[2]
op_file_ccf_regions = sys.argv[3]

region_names_dict = {}

# get number of regions

for pid in range(1,42,2):

    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    nis_id_str = str(apid).zfill(3)
    labels_csv_file = f'{ip_folder_labels}/allen_anno_data_{nis_id_str}.csv'
    dprint(labels_csv_file)
    region_names = []
    out_tissue = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])
            out_tissue.append(row[8])

    # record region names for beads that are located inside tissue region
    for idx, val in enumerate(out_tissue):
        if (val=='F'):
            region_names_dict[region_names[idx]] = 0

dprint("Number of regions identified: ", len(region_names_dict))

sorted_region_names_dict = OrderedDict(sorted(region_names_dict.items()))

for idx, name in enumerate(sorted_region_names_dict):
    print(idx, name)
    sorted_region_names_dict[name] = idx

with open(op_file_ccf_regions, 'w') as outfile:
    json.dump(sorted_region_names_dict, outfile)


