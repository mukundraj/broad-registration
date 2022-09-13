"""
Generates (optionally) aggregated region labels for regionally aggregated and normalized gene
expression vis data generation in s9d_region_agg.py. Generates labels for all
CCF region as opposed to s9c_labelgen.py, which does a few selected regions.
Should be run before s9d_region_agg.py

2022-07-12 - update - now also performs aggregation to (temporarily) bypass anndata install issue on gcp instance

Usage:

python s9d_labelgen.py
    inp: data root
    inp: path to folder with allen label ids
    inp: path to integrated bead coords + gene expression file
    i/o: path to interim folder for R script to write to 
    out: path to folder to write out region labels

Usage example:

python src/python/scripts/analysis/s9d_labelgen.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s8_raw_data/integrated_mats \
    /data_v3_nissl_post_qc/s9_analysis/s9d/interim \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels

Supplementary:

scp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/interim [IP]:~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/

Created by Mukund on 2022-07-07
"""

import sys
from produtils import dprint
from pathlib import Path
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.allen as allen
import csv
import subprocess

data_root = sys.argv[1]
label_data_folder = data_root+sys.argv[2]
ip_folder_counts = data_root+sys.argv[3]
io_folder_interim = data_root+sys.argv[4]
io_folder_labels = sys.argv[1]+sys.argv[5] # first as output then as input


# get number of regions
for pid in range(1, 208, 2):

    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment


    # get region names/labels
    nis_id_str = str(apid).zfill(3)
    labels_csv_file = f'{label_data_folder}/allen_anno_data_{nis_id_str}.csv'
    dprint(labels_csv_file)
    region_ids = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row[8]=='F':
                region_ids.append(int(row[9]))
            elif row[8]=='T':
                region_ids.append(0)
            else:
                assert(False)



    dprint(len(region_ids))
    aggr_labels_csv_file = f'{io_folder_labels}/agg_labels_{nis_id_str}.csv'
    with open(aggr_labels_csv_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, row in enumerate(region_ids):
            line = [row]
            writer.writerow(line)


    # perform aggregation in R
    subprocess.run(["Rscript", "src/python/scripts/analysis/s9d_prep.R", \
                    io_folder_labels, ip_folder_counts, io_folder_interim, \
                    str(apid)])
