"""
Script to generate hierarchical json data for dendrogram component in viewer.

Usage:

python  s1a_data_test.py \
    inp: data root

Usage example:

python src/python/scripts/analysis_sc/s1a_data_test.py \
    ~/Desktop/work/data/mouse_atlas \

Created by Mukund on 2022-08-11

"""

from produtils import dprint
import sys
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix

data_root = sys.argv[1]


ip_data_file1 = data_root+'/single_cell/s0_raw/jlanglie_scp_17_clusters_NZcounts.h5ad'
ip_data_file2 = data_root+'/single_cell/s0_raw/jlanglie_scp_newallAnnotAvgs.h5ad'
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


dprint("done")





