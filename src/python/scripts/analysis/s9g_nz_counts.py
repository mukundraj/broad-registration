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

python src/python/scripts/analysis/s9g_nz_counts.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels \
    /data_v3_nissl_post_qc/s8_raw_data/integrated_mats \
    /data_v3_nissl_post_qc/s9_analysis/s9g/nz_aggr_counts \

python src/python/scripts/analysis/s9g_nz_counts.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels_cshl \
    /data_v3_nissl_post_qc/s8_raw_data/integrated_mats \
    /data_v3_nissl_post_qc/s9_analysis/s9g/nz_aggr_counts_cshl_230222 \

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

pids = list(range(1, 208, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

for pid in pids:
# for pid in range(1, 4, 2):
    assert(pid!=5 and pid!=77 and pid!=167)

    # perform aggregation in R
    subprocess.run(["Rscript", "src/python/scripts/analysis/s9g_prep.R", \
                    ip_folder_labels, ip_folder_counts, out_folder, \
                    str(pid)])
