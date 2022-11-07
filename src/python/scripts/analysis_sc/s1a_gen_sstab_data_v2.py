
"""
Generates zarr files (nonzero counts and avg expr) for SingleCell viewer tab from version of data (csv files) generated on 2022-09-23. Also including celltype metadata on 2022-11-07.

Usage:

python  s1a_data_test.py
    inp: data root
    inp: path to nonzero counts file
    inp: path to avg vals file
    inp: path to cluster metadata file (for numcells for each cluster)
    inp: path to celltype metadata file
    out: output path

Usage example:

python src/python/scripts/analysis_sc/s1a_gen_sstab_data_v2.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_nz_mtx.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_avg_mtx.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/clusterSize.csv \
    /single_cell/s0/raw_v2/snRNA-seq_metadata.csv \
    /single_cell/s1 \

Supplementary:

// gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/single_cell/s0/zarr/scZarr.zarr gs://ml_portal/test_data
gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr.zarr gs://bcdportaldata/singlecell_data/scZarr.zarr

Created by Mukund on 2022-09-27

"""
from produtils import dprint
import sys
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix
import zarr
import pickle
import csv
import numcodecs
import re

data_root = sys.argv[1]
nz_csv_file = data_root+sys.argv[2]
avg_csv_file = data_root+sys.argv[3]
clustersize_csv_file = data_root+sys.argv[4]
metadata_file = data_root+sys.argv[5]
op_path = sys.argv[6]


# read metadata file
metadata = {}
with open(metadata_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    next(reader) # skip header
    for row in reader:
        # dprint(row[0], row[2], row[7], row[8])
        metadata[row[0]] = [row[2], row[7], row[8]] # [celltype, top_level_region, pct of cells in cluster from tlr]

dprint('metadata length:', len(metadata))




nClusters = 5030
nGenes = 21899

# nGenesTmp = 5

# create new zarr file
zarr_file = f'{data_root}/{op_path}/scZarr.zarr'
store = zarr.DirectoryStore(zarr_file) # https://zarr.readthedocs.io/en/stable/tutorial.html#storage-alternatives
root = zarr.group(store=store, overwrite=True)
nz_group = root.create_group(f'nz', overwrite=True)
nz_groupX = nz_group.zeros('X', shape=(nClusters, nGenes), chunks=(nClusters, 1), dtype='f4')
# nz_groupX[:] = 1

avg_group = root.create_group(f'avg', overwrite=True)
avg_groupX = avg_group.zeros('X', shape=(nClusters, nGenes), chunks=(nClusters, 1), dtype='f4')


obs_group = root.create_group(f'obs', overwrite=True)
var_group = root.create_group(f'var', overwrite=True)

metadataGroup = root.create_group(f'metadata', overwrite=True)

clusterSizes = []
clusterNames = []

# first, read clusterSize file
with open(clustersize_csv_file, 'r') as f:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(f)
    # Iterate over each row in the csv using reader object
    for idx,row in enumerate(csv_reader):
        # row variable is a list that represents a row in csv
        if (idx>0):
            clusterSizes.append(int(row[1]))
            clusterNames.append(row[0])



    # clustersArray = obs_group.zeros('clusters', shape=(nClusters), dtype='object', object_codec=numcodecs.JSON())
    clustersArray = obs_group.zeros('clusters', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
    clustersArray[:] = clusterNames


# next, create metadata arrays
cell_classes = [None]*nClusters
top_regions = [None]*nClusters
max_pcts = [None]*nClusters
# dprint(top_regions)

for idx, cname in enumerate(clusterNames):
    cname =  cname.split('=')[1]
    if cname in metadata:
        cell_classes[idx] = metadata[cname][0]
        top_regions[idx] = metadata[cname][1]
        max_pcts[idx] = metadata[cname][2][:-1] # remove % sign
    else:
        cell_classes[idx] = 'NA'
        top_regions[idx] = 'NA'
        max_pcts[idx] =  '0.0'

cellClassesArray = metadataGroup.zeros('cellclasses', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
cellClassesArray[:] = cell_classes
topRegionsArray = metadataGroup.zeros('topregions', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
topRegionsArray[:] = top_regions
maxPctsArray = metadataGroup.zeros('maxpcts', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
maxPctsArray[:] = max_pcts

uniqCellClasses = list(set(cell_classes))
uniqCellClassesArray = metadataGroup.zeros('uniqcellclasses', shape=(len(uniqCellClasses)), dtype='object', object_codec=numcodecs.VLenUTF8())
uniqCellClassesArray[:] = uniqCellClasses


nz_mat = np.zeros(shape=(nClusters, nGenes))
# next, populate nz data in zarr file
with open(nz_csv_file, 'r') as f:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(f)
    # Iterate over each row in the csv using reader object
    for idx,row in enumerate(csv_reader):
        # row variable is a list that represents a row in csv
        if (idx>0):
            curcell_nz_counts = np.array([int(x) for x in row[1:]])
            dprint('nzcsv', idx)
            nohead_idx = idx-1
            # nz_groupX[nohead_idx, :] = curcell_nz_counts[:nGenes] / clusterSizes[nohead_idx]
            nz_mat[nohead_idx, :] = curcell_nz_counts[:nGenes] / clusterSizes[nohead_idx]


        if idx==0:
            geneNames = row[1:]
            geneNamesArray = var_group.zeros('genes', shape=(nGenes), dtype='object', object_codec=numcodecs.VLenUTF8())
            # pgroup_genes_arr[0, :] = np.asarray(gene_names)
            geneNamesArray[:] = np.asarray(geneNames)
            p = re.compile("(.+?)=.+$")
            geneNamesArray[:] = np.asarray([p.search(x).group(1) for x in geneNamesArray])

nz_groupX[:] = nz_mat

avg_mat = np.zeros(shape=(nClusters, nGenes))

# third, populate avgs data in zarr file
with open(avg_csv_file, 'r') as f:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(f)
    # Iterate over each row in the csv using reader object
    for idx,row in enumerate(csv_reader):
        # row variable is a list that represents a row in csv
        if (idx>0):
            curcell_avg_counts = np.array([float(x) for x in row[1:]])
            dprint('avgs', idx)
            nohead_idx = idx-1
            avg_mat[nohead_idx, :] = curcell_avg_counts[:nGenes]

avg_groupX[:] = avg_mat


z = zarr.open(zarr_file)
dprint(z.tree())
# z.avg.X[:5, :2] z.nz.X[:5, :2] exit(0)
