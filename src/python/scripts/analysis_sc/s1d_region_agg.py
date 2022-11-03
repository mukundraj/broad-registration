"""
Generate regionally aggregated cell score data using R script ./s1d_prep.R

Usage:

python s1d_region_agg.py \
    inp: data root
    inp: path to processed region labels created by s9d_labelgen.py
    inp: path to scores data in anndata/h5ad format
    inp inp: start_pid and end_pid
    out: path to output folder for R script to write into

Usage example:

python src/python/scripts/analysis_sc/s1d_region_agg.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels \
    /single_cell/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix \
    1 207 \
    /single_cell/s1/s1d_region_agg \

Created by Mukund on 2022-10-28

"""

import sys
from produtils import dprint
import subprocess
import time
from multiprocessing import Pool

data_root = sys.argv[1]
ip_folder_labels = data_root+sys.argv[2]
ip_folder_scores = data_root+sys.argv[3]
start_pid = sys.argv[4]
end_pid = sys.argv[5]
op_folder = data_root+sys.argv[6]


def process_pid(pid):

    dprint(f'starting pid {pid}..................')
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # perform aggregation in R
    subprocess.run(["Rscript", "src/python/scripts/analysis_sc/s1d_prep.R", \
                    ip_folder_labels, ip_folder_scores, op_folder, \
                    str(apid)])

pids = list(range(int(start_pid), int(end_pid)+1, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)
if __name__ == '__main__':
    start = time.time()
    with Pool(1) as p:
        p.map(process_pid, pids)
    end = time.time()
    dprint(f'Total time {end - start} seconds.')
