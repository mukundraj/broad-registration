"""
Generates zarr files (nonzero counts and avg expr) for SingleCell viewer tab from version of data (csv files) generated on 2022-09-23. Also including celltype metadata on 2022-11-07.

Usage:

python  s1a_data_test.py
    inp: data root
    inp: path to nonzero counts file
    inp: path to avg vals file
    inp: path to cluster metadata file (for numcells for each cluster)
    inp: path to celltype metadata file
    inp: path to metadata with celltype cluster [added on 2022-12-08, updated on 2023-02-06]
    inp: path to CellSpatial tab's score histogram data
    inp: path to sum counts matrix  [avg is this val / numcells in cluster]
    out: output path

Usage example:

python src/python/scripts/analysis_sc/s1a_gen_sstab_data_v2.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_nz_mtx.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_avg_mtx.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/clusterSize.csv \
    /single_cell/s0/raw_v2/snRNA-seq_metadata.csv \
    /single_cell/s0/raw_v2/celltype_metadata/data_MouseAtlas_Submission_CellType_Metadata.tsv \
    /cell_spatial/s2/s2c/cell_jsons_s2c \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_sumCounts_mtx.csv \
    /single_cell/s0/raw_v2/neuropeptide_data \
    /single_cell/s1/scZarr_230221.zarr \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_230221.zarr gs://bcdportaldata/batch_230131/singlecell_data/scZarr_230221.zarr

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_230221.zarr/metadata gs://bcdportaldata/batch_230131/singlecell_data/scZarr_230221.zarr/metadata

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
import os

data_root = sys.argv[1]
nz_csv_file = data_root+sys.argv[2]
avg_csv_file = data_root+sys.argv[3]
clustersize_csv_file = data_root+sys.argv[4]
metadata_file = data_root+sys.argv[5]
celltype_metadata_file = data_root+sys.argv[6]
hist_data_path = data_root+sys.argv[7]
counts_csv_file = data_root+sys.argv[8]
neuropeptide_data_path = data_root+sys.argv[9]
op_zarr = sys.argv[10]


# read metadata file
metadata = {}
with open(metadata_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    next(reader) # skip header
    for row in reader:
        # dprint(row[0], row[2], row[7], row[8])
        metadata[row[0]] = [row[2], row[7], row[8]] # [celltype, top_level_region, pct of cells in cluster from tlr]

dprint('metadata length:', len(metadata))

# read top structure metadata file and populate to dict mapping cellname to structure
ctype_to_struct = {}
neuroTrans = {}
with open(celltype_metadata_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader) # skip header
    Max_TopStruct_idx = header.index('Max_TopStruct') # rank 1
    # Imputed_Max_TopStruct_idx = header.index('Imputed_Max_TopStruct')
    # Mapped_Max_TopStruct = header.index('Mapped_Max_TopStruct') 
    NumBeadsConfMapped_idx = header.index('NumBeadsConfMapped') # tells us if imputed or not - if nonzero number of beads confidently mapped
    NT_binary_idx = header.index('NT_binary') # index of neuropep binary
    for row in reader:
        # dprint(row)
        # ctype_to_struct[row[1]] = row[0]

        ctype_to_struct[row[0]] = '-'
        if (row[NumBeadsConfMapped_idx] != ''):
            ctype_to_struct[row[0]] = row[Max_TopStruct_idx]
        else:
            ctype_to_struct[row[0]] = row[Max_TopStruct_idx] + ' (imputed)'
            if (len(row[Max_TopStruct_idx])) == 0:
                ctype_to_struct[row[0]] = '-' # if imputed and yet no structure, then set to '-'

        if (row[NT_binary_idx] != ''): # if neuroTrans is not empty
            neuroTrans[row[0]] = row[NT_binary_idx]
        else:
            neuroTrans[row[0]] = '-'

dprint('struct metadata length:', len(ctype_to_struct))
# dprint('struct metadata:', ctype_to_struct)

# print ctype_to_struct where values in list are not equal
# for k, v in ctype_to_struct.items():
#     if v == '-':
#         dprint(k, v)

neuroPep = {}
neuroPepReceps = {}
# read in neuro peptide csv data
with open(neuropeptide_data_path+'/data_MouseAtlas_Submission_COMBINED_np_making_thresh.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    # header = next(reader) # skip header
    for row in reader:
        if (row[1] != ''):
            neuroPep[row[0]] = row[1]
        else:
            neuroPep[row[0]] = '-'

# read neuro peptide receptor csv data
with open(neuropeptide_data_path+'/data_MouseAtlas_Submission_COMBINED_np_sensing_thresh.tsv', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    # header = next(reader) # skip header
    for row in reader:
        if (row[1] != ''):
            neuroPepReceps[row[0]] = row[1]
        else:
            neuroPepReceps[row[0]] = '-'


nClusters = 5030
nGenes = 21899

# nGenesTmp = 5

# create new zarr file
zarr_file = f'{data_root}/{op_zarr}'
store = zarr.DirectoryStore(zarr_file) # https://zarr.readthedocs.io/en/stable/tutorial.html#storage-alternatives
root = zarr.group(store=store, overwrite=True)

nz_pct_group = root.create_group(f'nz_pct', overwrite=True)
nz_pct_groupX = nz_pct_group.zeros('X', shape=(nClusters, nGenes), chunks=(nClusters, 1), dtype='f4')
# nz_groupX[:] = 1

# create nx_counts_group
counts_group = root.create_group('counts', overwrite=True)
counts_groupX = counts_group.zeros('X', shape=(nClusters, nGenes), chunks=(nClusters, 1), dtype='f4') 

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
top_structs = [None]*nClusters
max_pcts = [None]*nClusters
map_status = [None]*nClusters # whether the cluster is mapped to a spatial celltype
neuropep = [None]*nClusters # neuropeptide
neuropep_recep = [None]*nClusters # neuropeptide receptor
neurotrans = [None]*nClusters # neurotransmitter binary

for idx, cname in enumerate(clusterNames):
    cname =  cname.split('=')[1]
    mapFilename = f'{hist_data_path}/{cname}.json'
    cnameMapExists = 'Y' if os.path.isfile(mapFilename)==True else 'N'
    if cname in metadata:
        cell_classes[idx] = metadata[cname][0]
        # top_regions[idx] = metadata[cname][1]
        max_pcts[idx] = metadata[cname][2][:-1] # remove % sign
        map_status[idx] = cnameMapExists
    else:
        cell_classes[idx] = 'NA'
        # top_regions[idx] = 'NA'
        max_pcts[idx] =  '0.0'
        map_status[idx] = cnameMapExists

    if cname in ctype_to_struct:
        top_structs[idx] = ctype_to_struct[cname]
    else:
        top_structs[idx] = '-'

    if cname in neuroTrans:
        neurotrans[idx] = neuroTrans[cname]
        top_structs[idx] += f' | {neuroTrans[cname]}'
    else:
        neurotrans[idx] = '-'

    if cname in neuroPep:
        neuropep[idx] = neuroPep[cname]
        top_structs[idx] += f' | {neuroPep[cname]}'
    else:
        neuropep[idx] = '-'

    if cname in neuroPepReceps:
        neuropep_recep[idx] = neuroPepReceps[cname]
        top_structs[idx] += f' | {neuroPepReceps[cname]}'
    else:
        neuropep_recep[idx] = '-'


cellClassesArray = metadataGroup.zeros('cellclasses', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
cellClassesArray[:] = cell_classes
topStructsArray = metadataGroup.zeros('topstructs1', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
topStructsArray[:] = top_structs
maxPctsArray = metadataGroup.zeros('maxpcts', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
maxPctsArray[:] = max_pcts
mapStatusArray = metadataGroup.zeros('mapStatus', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
mapStatusArray[:] = map_status

uniqCellClasses = list(set(cell_classes))
uniqCellClassesArray = metadataGroup.zeros('uniqcellclasses', shape=(len(uniqCellClasses)), dtype='object', object_codec=numcodecs.VLenUTF8())
uniqCellClassesArray[:] = uniqCellClasses

neuroPepArray = metadataGroup.zeros('neuroPep', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
neuroPepArray[:] = neuropep
neuroPepRecepArray = metadataGroup.zeros('neuroPepRecep', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
neuroPepRecepArray[:] = neuropep_recep
neuroTransArray = metadataGroup.zeros('neuroTrans', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
neuroTransArray[:] = neurotrans

nz_pct_mat = np.zeros(shape=(nClusters, nGenes))


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
            nz_pct_mat[nohead_idx, :] = curcell_nz_counts[:nGenes] / clusterSizes[nohead_idx]

        if idx==0:
            geneNames = row[1:]
            geneNamesArray = var_group.zeros('genes', shape=(nGenes), dtype='object', object_codec=numcodecs.VLenUTF8())
            # pgroup_genes_arr[0, :] = np.asarray(gene_names)
            geneNamesArray[:] = np.asarray(geneNames)
            p = re.compile("(.+?)=.+$")
            geneNamesArray[:] = np.asarray([p.search(x).group(1) for x in geneNamesArray])

nz_pct_groupX[:] = nz_pct_mat

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

counts_mat = np.zeros(shape=(nClusters, nGenes))
# fourth, populate counts csv
with open(counts_csv_file, 'r') as f:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(f)
    # Iterate over each row in the csv using reader object
    for idx,row in enumerate(csv_reader):
        # row variable is a list that represents a row in csv
        if (idx>0):
            curcell_counts = np.array([int(x) for x in row[1:]])
            dprint('counts', idx)
            nohead_idx = idx-1
            counts_mat[nohead_idx, :] = curcell_counts[:nGenes]


counts_groupX[:] = counts_mat

globalMaxAvgVal = str(round(np.max(avg_mat)))
dprint('globalMaxAvgVal', globalMaxAvgVal)
globalMaxAvgValArray = metadataGroup.zeros('globalMaxAvgVal', shape=(1), dtype='object', object_codec=numcodecs.VLenUTF8())
globalMaxAvgValArray[:] = globalMaxAvgVal

dprint('zarr file', zarr_file)
z = zarr.open(zarr_file)
dprint(z.tree())
# z.avg.X[:5, :2] z.nz.X[:5, :2] exit(0)
