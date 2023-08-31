"""Script to generate cell metadata file for NIH

Usage:
    io: data root
    ip: files containing requested genes list from NIH
    ip: path to zarr file with metadata
    op: cell metadata file

Example:

python src/python/scripts/misc/for_nih/get_cell_metadata.py \
    ~/Desktop/work/data/mouse_atlas \
    /misc/for_nih/s0/mouse_GPCR_Targets.csv \
    /single_cell/s1/scZarr_230321.zarr \
    /misc/for_nih/s1/celltype_metadata.csv \

Created by Mukund on 2023-08-31
"""

import zarr
import csv
import sys
from produtils import dprint
import numpy as np
import re

data_root = sys.argv[1]
gene_list_file = data_root+sys.argv[2]
metadata_zarr_file = data_root+sys.argv[3]
out_csv_file = data_root+sys.argv[4]

# get genes list
genes_list = []
# read in gene list csv file and ignore first and second rows
with open(gene_list_file, 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        genes_list.append(row[3])

genes_list = genes_list[2:]

geneNames = None
# read avg csv file
geneNames = None


z=zarr.open(metadata_zarr_file, mode='r')

cellclasses = z.metadata.cellclasses[:]
genecovers = z.metadata.genecovers[:]
neuropep = z.metadata.neuroPep[:]
neuropeprecep = z.metadata.neuroPepRecep[:]
neurotrans = z.metadata.neuroTrans[:]
topstructs = z.metadata.topstructs1[:]

clusterNames = z.obs.clusters[:]

# preserve only part after = in clusterNames
for idx,clusterName in enumerate(clusterNames):
    clusterNames[idx] = clusterName.split('=')[1]



# write all cell metadata to csv file
with open(out_csv_file, 'w') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(['SrNo', 'CellCluster','CellClass', 'BrainRegion', 'DefiningGeneSet','NeuroTrans','NeuroPep', 'NeuroPepRecep' ])
    for idx in range(len(cellclasses)):
        writer.writerow([idx+1, clusterNames[idx], cellclasses[idx], topstructs[idx], genecovers[idx], neurotrans[idx], neuropep[idx], neuropeprecep[idx]])
