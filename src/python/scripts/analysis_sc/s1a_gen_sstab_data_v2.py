"""
Generates zarr files (nonzero counts and avg expr) for SingleCell viewer tab from version of data (csv files) generated on 2022-09-23. Also including celltype metadata on 2022-11-07. Added clade info on 2023-10-04.

Usage:

python  s1a_data_test.py
    ip: data root
    ip: path to nonzero counts file
    ip: path to avg vals file
    ip: path to cluster metadata file (for numcells for each cluster)
    ip: path to celltype metadata file
    ip: path to metadata with celltype cluster [added on 2022-12-08, updated on 2023-02-06]
    ip: path to CellSpatial tab's score histogram data
    ip: path to sum counts matrix  [avg is this val / numcells in cluster]
    ip: processed clade file
    ip: additional cluster metadata file
    op: output path

Usage example:

python src/python/scripts/analysis_sc/s1a_gen_sstab_data_v2.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_nz_mtx.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_avg_mtx.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/clusterSize.csv \
    /single_cell/s0/raw_v2/snRNA-seq_metadata.csv \
    /single_cell/s0/raw_v2/celltype_metadata/CellType_Metadata__withSetCover.tsv \
    /cell_spatial/s2/s2c/cell_jsons_s2c \
    /single_cell/s0/raw_v2/20220912_QC_summary/cluster_sumCounts_mtx.csv \
    /single_cell/s0/raw_v2/neuropeptide_data \
    /single_cell/s1/processed_clade_info_v2.csv \
    /single_cell/s0/additional_metadata.csv \
    /single_cell/s1/scZarr_240318.zarr \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr gs://bcdportaldata/batch_231112/single_cell/s1/scZarr_240318.zarr

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr/metadata gs://bcdportaldata/batch_231112/singlecell_data/scZarr_240318.zarr/metadata

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

from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import src.python.utils.io as io



data_root = sys.argv[1]
nz_csv_file = data_root+sys.argv[2]
avg_csv_file = data_root+sys.argv[3]
clustersize_csv_file = data_root+sys.argv[4]
metadata_file = data_root+sys.argv[5]
celltype_metadata_file = data_root+sys.argv[6]
hist_data_path = data_root+sys.argv[7]
counts_csv_file = data_root+sys.argv[8]
neuropeptide_data_path = data_root+sys.argv[9]
proc_clade_file = data_root+sys.argv[10]
addtl_metadata_file = data_root+sys.argv[11]
op_zarr = sys.argv[12]

additional_metadata_dict = io.get_additional_cluster_metadata(addtl_metadata_file)

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
geneCovers = {}
with open(celltype_metadata_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader) # skip header
    Max_TopStruct_idx = header.index('Max_TopStruct') # rank 1
    # Imputed_Max_TopStruct_idx = header.index('Imputed_Max_TopStruct')
    # Mapped_Max_TopStruct = header.index('Mapped_Max_TopStruct') 
    NumBeadsConfMapped_idx = header.index('NumBeadsConfMapped') # tells us if imputed or not - if nonzero number of beads confidently mapped
    NT_binary_idx = header.index('NT_binary') # index of neuropep binary

    GeneCover_idx = header.index('Displayed_Gene_List_Cover') # index of gene list cover

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

        if (row[GeneCover_idx] != ''): # if gene cover is not empty
            geneCovers[row[0]] = row[GeneCover_idx].replace('|', '_') # to avoid mix up with | in SC table column
        else:
            geneCovers[row[0]] = '-'

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

portalClusterNames = []
indicesToRemove = []

clade_names = []
clade_annotations = []
# read proc_clade_file
with open(proc_clade_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if (row[1]=='Clade0' or row[2]==''):
            clade_names.append('-')
            indicesToRemove.append(len(clade_names)-1)
        else:
            clade_names.append(row[1]) # check for Clade0
        # clade_names.append(row[1]) # check for Clade0
        if (row[1]=='MC_37' or row[1]=='MC_39'):
            # Change clade names based on Jonah S's request on 230313
            clade_annotations.append('Isocortical Neurons')
        else:
            clade_annotations.append(row[2])

        # if (row[1]=='Clade0' or row[2]==''): # remove Clade0 and clades with no human readable name
        #     indicesToRemove.append(len(clade_names)-1)


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




# create a new list of cluster names (portalClusterNames) removing 
# - cells missing in metadata (meadata_file)
# - cells present in metadata but with cellclass = ''
# - metadata with no human readable clade name (clade_annotation)
# also set the  nClusters as len of portalClusterNames
# also make sure cladesArray and cellclasses array are of same length as portalClusterNames


# first, remove cells missing in metadata
for cluster in clusterNames: 
    cname =  cluster.split('=')[1]
    if cname in metadata:
        if metadata[cname][0] == '':
            indicesToRemove.append(clusterNames.index(cluster))
    else:
        indicesToRemove.append(clusterNames.index(cluster))

dprint('indicesToRemove', len(indicesToRemove))
dprint('len unique indicesToRemove', len(set(indicesToRemove)))
uniqIndicesToRemove = list(set(indicesToRemove))

# portalClusterNames by removing indicesToRemove from clusterNames
portalClusterNames = [i for j, i in enumerate(clusterNames) if j not in uniqIndicesToRemove]
portalCladeAnnotations = [i for j, i in enumerate(clade_annotations) if j not in uniqIndicesToRemove]
portalCladeNames = [i for j, i in enumerate(clade_names) if j not in uniqIndicesToRemove]

# nClusters = 5030
nClusters = len(clusterNames) - len(uniqIndicesToRemove)
nGenes = 21899

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

cladesArray = metadataGroup.zeros('clades', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
cladesArray[:] = portalCladeNames

cladeAnnotationsArray = metadataGroup.zeros('cladeAnnotations', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
cladeAnnotationsArray[:] = portalCladeAnnotations



# clustersArray = obs_group.zeros('clusters', shape=(nClusters), dtype='object', object_codec=numcodecs.JSON())
clustersArray = obs_group.zeros('clusters', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
clustersArray[:] = portalClusterNames


# next, create metadata arrays
cell_classes = [None]*nClusters
top_structs = [None]*nClusters
max_pcts = [None]*nClusters
map_status = [None]*nClusters # whether the cluster is mapped to a spatial celltype
neuropep = [None]*nClusters # neuropeptide
neuropep_recep = [None]*nClusters # neuropeptide receptor
neurotrans = [None]*nClusters # neurotransmitter binary
genecovers = [None]*nClusters # gene covers
additional_metadata = [None]*nClusters # additional metadata


for idx, cname in enumerate(portalClusterNames):
    cname =  cname.split('=')[1]
    mapFilename = f'{hist_data_path}/{cname}.json'
    cnameMapExists = 'Y' if os.path.isfile(mapFilename)==True else 'N'
    if cname in metadata:
        cell_classes[idx] = metadata[cname][0]
        # top_regions[idx] = metadata[cname][1]
        max_pcts[idx] = metadata[cname][2][:-1] # remove % sign
        map_status[idx] = cnameMapExists
    else:
        # cell_classes[idx] = 'NA'
        cell_classes[idx] = '-' # make empty if not found
        # top_regions[idx] = 'NA'
        max_pcts[idx] =  '0.0'
        map_status[idx] = cnameMapExists
        dprint(cname, 'not found in metadata file')

    if cname in ctype_to_struct:
        top_structs[idx] = ctype_to_struct[cname]
    else:
        top_structs[idx] = '-'

    if cname in geneCovers:
        # top_structs[idx] += f' | {geneCovers[cname]}'
        genecovers[idx] = geneCovers[cname]
    else:
        genecovers[idx] = '-'

    if cname in neuroTrans:
        neurotrans[idx] = neuroTrans[cname]
        # top_structs[idx] += f' | {neuroTrans[cname]}'
    else:
        neurotrans[idx] = '-'

    if cname in neuroPep:
        neuropep[idx] = neuroPep[cname]
        # top_structs[idx] += f' | {neuroPep[cname]}'
    else:
        neuropep[idx] = '-'

    if cname in neuroPepReceps:
        neuropep_recep[idx] = neuroPepReceps[cname]
        # top_structs[idx] += f' | {neuroPepReceps[cname]}'
    else:
        neuropep_recep[idx] = '-'

    if cname in additional_metadata_dict:
        additional_metadata[idx] = additional_metadata_dict[cname]
    else:
        additional_metadata[idx] = ''

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

geneCoversArray = metadataGroup.zeros('geneCovers', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
geneCoversArray[:] = genecovers

additionalMetadataArray = metadataGroup.zeros('additionalMetadata', shape=(nClusters), dtype='object', object_codec=numcodecs.VLenUTF8())
additionalMetadataArray[:] = additional_metadata



nClustersFull = 5030
nz_pct_mat = np.zeros(shape=(nClustersFull, nGenes))
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

# remove rows in nz_pct_mat that are in indicesToRemove
nz_pct_mat = np.delete(nz_pct_mat, uniqIndicesToRemove, axis=0)

nz_pct_groupX[:] = nz_pct_mat


avg_mat = np.zeros(shape=(nClustersFull, nGenes))

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

# remove rows in avg_mat that are in indicesToRemove
avg_mat = np.delete(avg_mat, uniqIndicesToRemove, axis=0)
avg_groupX[:] = avg_mat


counts_mat = np.zeros(shape=(nClustersFull, nGenes))
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


# remove rows in counts_mat that are in indicesToRemove
counts_mat = np.delete(counts_mat, uniqIndicesToRemove, axis=0)
counts_groupX[:] = counts_mat

globalMaxAvgVal = str(round(np.max(avg_mat)))
dprint('globalMaxAvgVal', globalMaxAvgVal)
globalMaxAvgValArray = metadataGroup.zeros('globalMaxAvgVal', shape=(1), dtype='object', object_codec=numcodecs.VLenUTF8())
globalMaxAvgValArray[:] = globalMaxAvgVal

dprint('zarr file', zarr_file)
z = zarr.open(zarr_file)
dprint(z.tree())
# z.avg.X[:5, :2] z.nz.X[:5, :2] exit(0)
