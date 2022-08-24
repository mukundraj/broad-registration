"""
Script to generate hierarchical json data for dendrogram component in viewer.

Usage:

python  s1a_data_test.py \
    inp: data root
    out: output path

Usage example:

python src/python/scripts/analysis_sc/s1a_data_test.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s0 \

Supplementary:

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/single_cell/s0/zarr/z1.zarr gs://ml_portal2/test_data2/single_cell/s0/zarr

Created by Mukund on 2022-08-11

"""

from produtils import dprint
import sys
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix
import zarr

data_root = sys.argv[1]
op_path = sys.argv[2]

ip_data_file1 = data_root+'/single_cell/s0/raw/jlanglie_scp_17_clusters_NZcounts.h5ad'
ip_data_file2 = data_root+'/single_cell/s0/raw/jlanglie_scp_newallAnnotAvgs.h5ad'
data1 = ann.read_h5ad(ip_data_file1)
data2 = ann.read_h5ad(ip_data_file2)

csr_data1 = csr_matrix(data1.X)
csr_data2 = csr_matrix(data2.X)

dprint('obs_names', data1.obs_names)
dprint('var_names', data1.var_names)

dprint('len obs_names', len(list(data1.obs_names)))
dprint('len var_names', len(list(data1.var_names)))
dprint("jlanglie_scp_17_clusters_NZpercent.h5ad maxmin", csr_data1.max(), csr_data1.min())
dprint("jlanglie_scp_newallAnnotAvgs.h5ad maxmin", csr_data2.max(), csr_data2.min())


zarr_file = f'{data_root}/{op_path}/zarr/z1.zarr'

data1.write_zarr(zarr_file, (4864,1)) # storing 1 gene per chunk

zfile = zarr.open(zarr_file, mode='r+') # no need to close explicitly - https://zarr.readthedocs.io/en/stable/tutorial.html

dprint(zarr_file)

dprint("done")

# z = zarr.open('/Users/{username}/Desktop/work/data/mouse_atlas//single_cell/s0/zarr/z1.zarr')
# z.tree()
# z.obs._index
# z.X[:5, :10]
# z.X.info to get chunk shape






