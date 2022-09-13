"""
Generates zarr files (nonzero counts and avg expr) for SingleCell viewer tab

Usage:

python  s1a_data_test.py
    inp: data root
    inp: path to nonzero counts file
    inp: path to avg vals file
    inp: path to cluster metadata file (for numcells for each cluster)
    out: output path

Usage example:

python src/python/scripts/analysis_sc/s1a_data_test.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s0/raw/jlanglie_scp_17_clusters_NZcounts.h5ad \
    /single_cell/s0/raw/jlanglie_scp_newallAnnotAvgs.h5ad \
    /single_cell/s0/raw/jlanglie_scp_18_cluster_meta_for_Mukund.pickle \
    /single_cell/s0 \

Supplementary:

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/single_cell/s0/zarr/z_proportions.zarr gs://ml_portal/test_data
gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/single_cell/s0/zarr/z_avgs.zarr gs://ml_portal/test_data

Created by Mukund on 2022-08-11

"""

from produtils import dprint
import sys
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix
import zarr
import pickle

data_root = sys.argv[1]
ip_data_file1 = data_root+sys.argv[2]
ip_data_file2 = data_root+sys.argv[3]
ip_data_metadata = data_root+sys.argv[4]
op_path = sys.argv[5]

with open(ip_data_metadata, 'rb') as f:
    x = pickle.load(f)
    dprint(x.columns)
    # dprint(x.iloc[0])
    # dprint(x.iloc[0,4])
    dprint(x["count"])
    counts = x["count"].values

dprint(counts)
data1 = ann.read_h5ad(ip_data_file1)
data2 = ann.read_h5ad(ip_data_file2)

data1.X = np.round(data1.X/counts[:,None],2)
data2.X = np.round(data2.X.astype(np.float64), 2)

csr_data1 = csr_matrix(data1.X)
csr_data2 = csr_matrix(data2.X)

dprint('obs_names', data1.obs_names)
dprint('var_names', data1.var_names)

dprint('len obs_names', len(list(data1.obs_names)))
dprint('len var_names', len(list(data1.var_names)))
dprint("jlanglie_scp_17_clusters_NZpercent.h5ad maxmin", csr_data1.max(), csr_data1.min())
dprint("jlanglie_scp_newallAnnotAvgs.h5ad maxmin", csr_data2.max(), csr_data2.min())

dprint(data2.X)
zarr_file_proportions = f'{data_root}/{op_path}/zarr/z_proportions.zarr'
data1.write_zarr(zarr_file_proportions, (4864,1)) # storing 1 gene per chunk
zfile_proportions = zarr.open(zarr_file_proportions, mode='r+') # no need to close explicitly - https://zarr.readthedocs.io/en/stable/tutorial.html


zarr_file_avgs = f'{data_root}/{op_path}/zarr/z_avgs.zarr'
data2.write_zarr(zarr_file_avgs, (4864,1)) # storing 1 gene per chunk
zfile_avgs = zarr.open(zarr_file_avgs, mode='r+') # no need to close explicitly - https://zarr.readthedocs.io/en/stable/tutorial.html

dprint("done")

# z = zarr.open('/Users/{username}/Desktop/work/data/mouse_atlas//single_cell/s0/zarr/z1.zarr')
# z.tree()
# z.obs._index
# z.X[:5, :10]
# z.X.info to get chunk shape






