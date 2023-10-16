"""
Query for cell cluster names for cluster classified as endo, fibroblast, and
micropahges

Usage:

phthon evan_query.py
    ip: data_root
    ip: scZarr path relative to data_root
    op: output folder

Example:

python tmp_scripts/evan_query.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s1/scZarr_231011.zarr \
    /misc/for_evan \

Created by Mukund on 2023-10-16
"""

import sys
import zarr
from produtils import dprint
import os


data_root = sys.argv[1]
zarr_relative = sys.argv[2]
out_folder = sys.argv[3]

zarr_path = data_root + zarr_relative
outfile = data_root + out_folder + '/output.csv'

dprint(zarr_path)

z = zarr.open(zarr_path, mode='r')

dprint(z.tree())
cellclasses = z.metadata['cellclasses']

# get indices of endo, fibroblast, and macrophage in cellclasses

endo_idx = [i for i, x in enumerate(cellclasses) if x == 'Endo']
fibro_idx = [i for i, x in enumerate(cellclasses) if x == 'fibro']
macr_idx = [i for i, x in enumerate(cellclasses) if x == 'macrophage']
dprint(len(endo_idx))
dprint(len(fibro_idx))
dprint(len(macr_idx))


clusters = z.obs['clusters']

# get clusters for indices in endo_idx, fibro_idx, and macr_idx


endo_clusters = [clusters[i] for i in endo_idx]
fibro_clusters = [clusters[i] for i in fibro_idx]
macr_clusters = [clusters[i] for i in macr_idx]


# keep only part after = in cluster names
endo_clusters = [x.split('=')[1] for x in endo_clusters]
fibro_clusters = [x.split('=')[1] for x in fibro_clusters]
macr_clusters = [x.split('=')[1] for x in macr_clusters]

# write out a csv with cluster name and cell class as columns for each cell class

with open(outfile, 'w') as f:
    f.write('cluster\tcellclass\n')
    for c in endo_clusters:
        f.write('endo\t'+c+'\n')
    for c in fibro_clusters:
        f.write('fibro\t'+c+'\n')
    for c in macr_clusters:
        f.write('macrophage\t'+c+'\n')








dprint(endo_clusters)


