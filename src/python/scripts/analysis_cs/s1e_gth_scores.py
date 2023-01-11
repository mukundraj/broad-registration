
"""
Generates regionwise aggregated, greater than threshod (assumes lt threshold is
zero) score matrix for cell scores. Also, generates a csv file with regionwise
total scores mapped to Allen regionids. Uses analysis_cs/s2e_prep.R for fast
aggregation R subroutines.

Usage:

python s1e_gth_scores.py
    inp: data root
    inp: path to folder with allen label ids
    inp: path to integrated bead coords + cell scores file
    out: path to output folder with nonzero cell scores

Usage example:

python src/python/scripts/analysis_cs/s1e_gth_scores.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels \
    /cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW \
    /cell_spatial/s1/s1e_gth_aggr_scores \

python src/python/scripts/analysis_cs/s1e_gth_scores.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels_cshl \
    /cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW \
    /cell_spatial/s1/s1e_gth_aggr_scores \

Created by Mukund on 2022-12-04

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
ip_folder_scores = data_root+sys.argv[3]
out_folder = data_root+sys.argv[4]

pids = list(range(1, 208, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

# pids = [1, 15]
for pid in pids:
# for pid in range(1, 4, 2):
    assert(pid!=5 and pid!=77 and pid!=167)

    # perform aggregation in R
    subprocess.run(["Rscript", "src/python/scripts/analysis_cs/s1e_prep.R", \
                    ip_folder_labels, ip_folder_scores, out_folder, \
                    str(pid)])

