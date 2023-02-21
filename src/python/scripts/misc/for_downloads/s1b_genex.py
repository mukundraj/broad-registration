"""Generate download files for gene expression tab data.

Usage:

python s1b_genex.py \
    ip: data root
    ip: gene expression matrix data - genes x beads
    ip: ccf labels for beads
    ip: bead coords in nissl space
    op: op folder

Example:

python src/python/scripts/misc/for_downloads/s1b_genex.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s8_raw_data/integrated_mats \
    /v3/s1/CSHL_CTMapping \
    /data_v3_nissl_post_qc/s4_bead_to_segid/bead_to_segid \
    /misc/for_downloads/s1b_genexp \


Created by Mukund on 2023-02-20
"""



from produtils import dprint
import anndata as ann
import sys
import os
import shutil
import zarr
import numpy as np
from scipy.sparse import csr_matrix



from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import src.python.utils.filename_maps as fmaps

data_root = sys.argv[1]
ip_integ_mats = f'{data_root}/{sys.argv[2]}'
ip_ccf_labels = f'{data_root}/{sys.argv[3]}'
ip_coords_nissp = f'{data_root}/{sys.argv[4]}'
op_folder = f'{data_root}/{sys.argv[5]}'

pidToSrno, srnoToPid = fmaps.get_pidToSrno()

# create a dicrectory scdata in opfolder and overwrite if exists
data_folder = f'{op_folder}/genex_data'
dprint(data_folder)
if os.path.exists(data_folder):
    shutil.rmtree(data_folder)
os.mkdir(data_folder)

# copy readmen.txt from s1d_single_cell folder to scdata folder
readme_file = f'{op_folder}/readme.txt'
shutil.copy(readme_file, data_folder)


for srno in range(1, 102):
# for srno in range(1, 3):
    pid = srnoToPid[str(srno).zfill(3)]

    ip_counts_file  = f'{ip_integ_mats}/ad_counts_{str(pid)}.h5ad'
    counts = ann.read_h5ad(ip_counts_file)
    # coords = ann.read_h5ad(ip_coords_file)
    ngenes, nbeads = counts.shape
    dprint(ngenes, nbeads)
    dprint(counts.X)

    an = ann.AnnData(X=counts.X)

    # add coords info to anndata file

    pid_str = str(pid).zfill(3)
    coords_file = f'{ip_coords_nissp}/bead_to_segid_{pid_str}.csv'

    coords = np.loadtxt(coords_file, delimiter=',', usecols=(2,3), dtype=np.int16)
    dprint(coords.shape)
    an.varm['coords'] = coords

    # add ccfRids to anndata file
    ccfRids = np.zeros(nbeads)

    ip_cc_labels_file = f'{ip_ccf_labels}/dd{pid}_CTMapping.h5ad'
    labels_file = ann.read_h5ad(ip_cc_labels_file)
    labels_inds = labels_file.var['rowName']
    labels_rids = labels_file.var['CSHL_CCFID']

    for i, ind in enumerate(labels_inds):
        int_ind = int(ind)-1 # R to python index
        ccfRids[int_ind] = labels_rids[i]
    an.var['ccfRids'] = ccfRids

    dprint(an)

    an_fname = f'{data_folder}/sr_{srno}.h5ad'
    # write out an_fname
    an.write(an_fname, compression="gzip")


