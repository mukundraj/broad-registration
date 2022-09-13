"""
Generates regionwise aggregated, nonzero count matrix. Also, generates a csv
file with regionwise total counts mapped to Allen regionids. Uses s2a_prep.R
for fast aggregation R subroutines.

Usage:

python s2a_nz_counts.py
    inp: data root
    inp: path to folder with allen label ids
    inp: path to integrated bead coords + gene expression file
    out: path to output folder with nonzero counts

Usage example:

python src/python/scripts/analysis_sc/s2a_nz_counts.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels \
    /data_v3_nissl_post_qc/s8_raw_data/integrated_mats \
    /single_cell/s2/nz_aggr_counts \

Created by Mukund on 2022-09-13
"""

import sys
from produtils import dprint
from pathlib import Path
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.allen as allen
import subprocess


data_root = sys.argv[1]
ip_folder_labels = data_root+sys.argv[2]
ip_folder_counts = data_root+sys.argv[3]
out_folder = data_root+sys.argv[4]


# for pid in range(1, 208, 2):
for pid in range(1, 4, 2):

    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # perform aggregation in R
    subprocess.run(["Rscript", "src/python/scripts/analysis_sc/s2a_prep.R", \
                    ip_folder_labels, ip_folder_counts, out_folder, \
                    str(apid)])
