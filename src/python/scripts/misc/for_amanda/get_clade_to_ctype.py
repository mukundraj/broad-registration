"""Generate a list of cell types in each clade when given a puckid.

python script.py \
    io: data root
    ip: srno
    ip: path to scZarr for celltype <-> clade mappings
    ip: path to aggr data folder to spatially determine celltypes <-> clades
    op: path to op csv file

python src/python/scripts/misc/for_amanda/get_clade_to_ctype.py \
    ~/Desktop/work/data/mouse_atlas \
    41 \
    /single_cell/s1/scZarr_240531.zarr \
    /cell_spatial/s1/s1c_aggr_data_240531 \
    /misc/for_amanda/clade_to_clusternames

Created by Mukund on 2025-04-17
"""
import sys

from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import src.python.utils.filename_maps as fmaps
import zarr
import anndata as an
import numpy as np

data_root = sys.argv[1]
srno = sys.argv[2]
scZarrFile = f'{data_root}{sys.argv[3]}'
ip_aggr_folder = f'{data_root}{sys.argv[4]}'


# pad srno with zero to have length 3
srno = str(srno).zfill(3)

pidToSrno, srnoToPid = fmaps.get_pidToSrno()

pid = srnoToPid[srno]
print(pid)

z = zarr.open(scZarrFile)
# print(z.metadata.clades)
# print(z.obs.clusters)
clades = list(z.metadata.clades)
clusters = list(z.obs.clusters)
clusters = [cluster.split('=')[1] for cluster in clusters]

cluster_to_clade = {}
for i, clade in enumerate(clades):
    cluster_to_clade[clusters[i]] = clade

print(cluster_to_clade)


padded_pid = str(pid).zfill(3)
print(ip_aggr_folder)

aggr_file = an.read_h5ad(f'{ip_aggr_folder}/dd{pid}_aggr_cladeclass.h5ad')
obs = aggr_file.obs['name'].tolist()
#
# get dims of aggr_file
m,n = aggr_file.shape
#
# # print aggr_file.obs without truncation
#
# how many celltypes, clades, and cellclasses are there?
celltypes_len = m - 240
celltypes = obs[:m-240]
# clades = obs[m-240:m-16]
# print(aggr_file.X.shape)
#
# # get submatrix of celltypes
celltypes_mat = aggr_file.X[:celltypes_len, :]
print(celltypes_mat.shape)

inds = celltypes_mat>0.3
r,c = inds.nonzero() # get row and column indices of nonzero elements


uniq_r = list(set(r))
# convert to int
uniq_r = [int(i) for i in uniq_r]
print(uniq_r)
print(len(uniq_r))

# clusternames from indices in uniq_r
clusternames = [clusters[i] for i in uniq_r]

cladenames = [cluster_to_clade[clusternames[i]] for i in range(len(clusternames))]
print(len(cladenames))
print(len(clusternames))

# create a clade to clusternames array map
clade_to_clusternames = {}
for i, clade in enumerate(cladenames):
    if clade not in clade_to_clusternames:
        clade_to_clusternames[clade] = []
    clade_to_clusternames[clade].append(clusternames[i])


# write out clade_to_clusternames in a output csv
op_csv = f'{data_root}{sys.argv[5]}/clade_to_clusternames_{srno}.csv'
with open(op_csv, 'w') as f:
    for clade in clade_to_clusternames:
        f.write(f'{clade},{",".join(clade_to_clusternames[clade])}\n')

