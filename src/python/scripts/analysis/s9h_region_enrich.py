"""
Region enrichment data generation. Generate zarr with allen region id to array
of num beads expressing gene (array length = number of genes).

Usage:

python s9h_region_enrich.py
    i/o: data root
    inp: nonzero, aggregated counts prepared by ./s2a_nz_counts.py
    out: path to output folder

Usage example:

python src/python/scripts/analysis/s9h_region_enrich.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9g/nz_aggr_counts \
    /data_v3_nissl_post_qc/s9_analysis/s9h \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9h/nz_aggr.zarr gs://bcdportaldata/cellspatial_data/nz_aggr.zarr

gsutil cp ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9h/gene_info.json gs://bcdportaldata/cellspatial_data/gene_info.json

Created by Mukund on 2022-09-12
"""


import sys
from produtils import dprint
from pathlib import Path
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.data_struct_helpers as dsh
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import nrrd
from produtils import dprint
import numpy as np
import anndata as ann
from scipy.sparse import csr_matrix, hstack
from treelib import Node, Tree
import zarr
import numcodecs
import json


data_root = sys.argv[1]
ip_nz_aggr_data_folder = data_root+sys.argv[2]
op_path = data_root+sys.argv[3]

reference_space_key = 'annotation/ccf_2017/'
# resolution = 10
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# rspc = ReferenceSpaceCache(resolution, reference_space_key)
# ID 1 is the adult mouse structure graph
allentree = rspc.get_structure_tree(structure_graph_id=1)
name_map = allentree.get_name_map()
ancestors_list = list(allentree.get_ancestor_id_map().values())

annotation, meta = rspc.get_annotation_volume()

# dprint(ancestors_list[0])

# io_folder_interim = "/Users/mraj/Desktop/work/data/mouse_atlas"+"/data_v3_nissl_post_qc/s9_analysis/s9d/interim"
# nis_id_str = "003"
# aggr_counts_file = f'{io_folder_interim}/aggr_counts_{nis_id_str}.h5ad'
# aggr_counts = ann.read_h5ad(aggr_counts_file)
# csr = csr_matrix(aggr_counts.X)
# dprint(csr)


# start

# tree = Tree()

# for nodes_list in ancestors_list:
#     dsh.add_nodes(tree, nodes_list, name_map)


# traverse tree

tree_nodes_info = []

# for hydrating tree nodes and populating an array of tree node data
def get_child_info(tree, child_id, name_map, parent_id, nGenes, rid_to_idx_local, dataXcsr, bead_counts):

    children = tree.children(child_id)
    data = {"rid":child_id, "nz_counts":np.zeros(nGenes), "num_beads":0}

    # if (child_id==1025 or child_id==10683):
    #     if child_id in rid_to_idx_local:
    #         # dprint('child_id', child_id, 'parent_id', parent_id)
    #         dprint('child_id', child_id, rid_to_idx_local)
    #         local_row_idx = rid_to_idx_local[child_id]
    #         dprint('local_row_idx', local_row_idx)
    #         dprint("child_id", child_id, 'nbeads',bead_counts[local_row_idx], 'counts', np.sum(np.squeeze(dataXcsr.getrow(local_row_idx).toarray())), 'numchildren', len(children))
    if (len(children)>0):
        if child_id in rid_to_idx_local:
            local_row_idx = rid_to_idx_local[child_id]
            if bead_counts[local_row_idx][1] > 0:
                data['nz_counts'] += np.squeeze(dataXcsr.getrow(local_row_idx).toarray())
                data['num_beads'] += bead_counts[local_row_idx][1]
        for cur_child_node in children:
            cur_child_id = cur_child_node.identifier
            info = get_child_info(tree, cur_child_id, name_map, child_id, nGenes, rid_to_idx_local, dataXcsr, bead_counts)

            # assimalate data into parent
            data['nz_counts'] += info['nz_counts']
            data['num_beads'] += info['num_beads']


    else: # leaf node - read data for region and update counts
        # data['nz_counts'] += np.zeros(3)
        if child_id in rid_to_idx_local:
            local_row_idx = rid_to_idx_local[child_id]
            data['nz_counts'] += np.squeeze(dataXcsr.getrow(local_row_idx).toarray())
            data['num_beads'] += bead_counts[local_row_idx][1]

    tree_nodes_info.append(data)
    return data

# for getting nGenes
ip_data_file1 = ip_nz_aggr_data_folder+"/nz_aggr_counts_003.h5ad"
data1 = ann.read_h5ad(ip_data_file1)
data1X = data1.X
data1Xcsr = csr_matrix(data1X)
# dprint(data1Xcsr)
# dprint(data1Xcsr.getrow(0))
_, nGenes = data1Xcsr.shape
# dprint("nGenes", nGenes)
# read total count file for mapping info from local index to rid

# also pass data1 and count file to get_child_info
# read puck data and get region ids
# nis_id_str = "003"
# aggr_counts_file = f'{ip_nz_aggr_data_folder}/nz_aggr_num_beads_{nis_id_str}.csv'
# bead_counts = np.genfromtxt(aggr_counts_file, delimiter=',', skip_header=1).astype(np.int32)

# rid_to_idx_local = {}
# for idx, x in enumerate(bead_counts):
#     rid_to_idx_local[x[0]] = idx

# m,n = np.shape(bead_counts)
# dprint('bead_counts shape :', m, n)

# data = get_child_info(tree, ancestors_list[0][0], name_map, "", nGenes, rid_to_idx_local, data1Xcsr)
# dprint(data)
# dprint(len(tree_nodes_info))
# dprint(tree_nodes_info[-3:])

# tree_nodes_info_idxed = [{ **x, 'idx':i} for i, x in enumerate(tree_nodes_info)]


zarr_filename = f'{op_path}/nz_aggr.zarr'
store = zarr.DirectoryStore(zarr_filename) # https://zarr.readthedocs.io/en/stable/tutorial.html#storage-alternatives
root = zarr.group(store=store, overwrite=True)


pids = list(range(1, 208, 2))
# pids = [1,207]
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

globalGroupX = None
globalGroupXout = None
gzVal = None
gzValOut = None

curPuckSums = None
maxExpr = None
maxExprPuck = None


for pid in pids:
    dprint("pid: ", pid)
    assert(pid!=5 and pid!=77 and pid!=167)
    # read regionwise gene exp data and regionwise total expression data
    nis_id_str = str(pid).zfill(3)
    ip_data_file = f'{ip_nz_aggr_data_folder}/nz_aggr_counts_{nis_id_str}.h5ad'
    # dprint(ip_data_file)
    data = ann.read_h5ad(ip_data_file)
    gene_names = data.var_names
    dataX = data.X
    dataXcsr = csr_matrix(dataX)
    M, N = dataXcsr.shape
    aggr_counts_file = f'{ip_nz_aggr_data_folder}/nz_aggr_num_beads_{nis_id_str}.csv'
    all_region_nz_array = np.squeeze(np.array(dataXcsr[1,:].sum(axis=0))) # excluding region zero (non tissue)

    bead_counts = np.genfromtxt(aggr_counts_file, delimiter=',', skip_header=1).astype(np.int32)
    m,n = np.shape(bead_counts)
    assert(m==M)
    rid_to_idx_local = {}
    for idx, x in enumerate(bead_counts):
        rid_to_idx_local[x[0]] = idx

    # create and hydrate tree based on current puck's data
    tree = Tree()
    for nodes_list in ancestors_list:
        dsh.add_nodes(tree, nodes_list, name_map)

    tree_nodes_info = []
    data = get_child_info(tree, ancestors_list[0][0], name_map, "", nGenes, rid_to_idx_local, dataXcsr, bead_counts)
    dprint("final tree_nodes_info ", tree_nodes_info[-1])
    rid_to_idx_map = {}
    rids = []
    for idx, x in enumerate(tree_nodes_info):
        rid_to_idx_map[x['rid']]=idx
        rids.append(x['rid'])
    nRids = len(tree_nodes_info)
    if pid==1:
        rids = np.array(rids)
        dprint(rids)
        rids_group = root.create_group(f'rids', overwrite=True)
        rids_groupX = rids_group.zeros('X', shape=(1,nRids), chunks=(1, nRids), dtype='i4')
        rids_groupX[0, :] = rids

        # for storing global (as opposed to puckwise) info
        globalGroup = root.create_group(f'pall', overwrite=True)
        globalGroupX = globalGroup.zeros('X', shape=(nRids, nGenes), chunks=(1, nGenes), dtype='f4')
        globalGroupX[:] = 0
        globalGroupXout = globalGroup.zeros('Xout', shape=(nRids, nGenes), chunks=(1, nGenes), dtype='f4')
        # globalGroupXout = np.tile(all_region_nz_array, (nRids, 1))  # initializing here to prevent multiple inclusion of exterior beads
        gzVal = np.zeros(nRids)
        gzValOut = np.zeros(nRids)

        maxExpr = np.zeros(nGenes)
        maxExprPuck = np.ones(nGenes)


    # write dict info to zarr
    # pgroup = root.create_group(f'p{pid}', overwrite=True)
    # pgroupX = pgroup.zeros('X', shape=(nRids, nGenes), chunks=(1, nGenes), dtype='f4')
    # pgroupX[:] = 0
    # pgroupXout = pgroup.zeros('Xout', shape=(nRids, nGenes), chunks=(1, nGenes), dtype='f4')
    # pgroupXout[:] = 0

    pgroupX = np.zeros((nRids, nGenes))
    pgroupXout = np.zeros((nRids, nGenes))

    all_region_beads = np.sum(np.array([x['num_beads'] for x in tree_nodes_info]), axis=0)
    # globalGroupXout += np.tile(all_region_nz_array_hydrated, (nRids, 1))  # initializing here to prevent multiple inclusion of exterior beads
    # dprint(np.shape(all_region_nz_array_hydrated))
    # exit(0)

    for rid in rid_to_idx_map: 
        global_region_idx = rid_to_idx_map[rid]

        zval = tree_nodes_info[global_region_idx]['num_beads']
        gzVal[global_region_idx] += zval

        # pgroupX[global_region_idx,:] = tree_nodes_info[global_region_idx]['nz_counts']
        globalGroupX[global_region_idx,:] += tree_nodes_info[global_region_idx]['nz_counts'] # will be normalized later
        pgroupX[global_region_idx,:] = tree_nodes_info[global_region_idx]['nz_counts']/zval # normalized her to determine puck with max expr per 10k beads

        # now calculate % nonzero outsize region
        zvalOut = tree_nodes_info[rid_to_idx_map[997]]['num_beads'] - tree_nodes_info[global_region_idx]['num_beads'] # denomminator -> num of beads outside current region
        gzValOut[global_region_idx] += zvalOut
        # if (gzValOut[global_region_idx] > 0):
        #     dprint(gzValOut[global_region_idx])

        # pgroupXout[global_region_idx,:] = all_region_nz_array - tree_nodes_info[global_region_idx]['nz_counts'] # numerator -> num of nonzero count beads outside current region
        pgroupXout[global_region_idx,:] = tree_nodes_info[rid_to_idx_map[997]]['nz_counts'] - tree_nodes_info[global_region_idx]['nz_counts'] # numerator -> num of nonzero count beads outside current region, will be normalized later
        globalGroupXout[global_region_idx,:] += pgroupXout[global_region_idx,:]
        pgroupXout[global_region_idx,:] /= zvalOut # normalized here to maintain consistency with pgroupX

    # dprint("DEBUG", pid, globalGroupX[1098, 8198], gzVal[1098], globalGroupXout[1098, 8198], gzValOut[1098])
    dprint("DEBUG", pid, tree_nodes_info[1098]['nz_counts'][8198] , gzVal[1098], pgroupXout[1098, 8198], gzValOut[1098])
    sys.stdout.flush()

    curPuckSums = np.sum(np.nan_to_num(pgroupX[:,:]), axis=0) # summmation over regions of regionwise normalized expression counts for all genes in current puck
    # dprint('pgroupX ', np.nan_to_num(pgroupX[:]))
    # dprint('curPuckSums', curPuckSums)
    # dprint(' max curPuckSums', np.max(curPuckSums))
    curGreaterIdxs = np.where(curPuckSums > maxExpr)
    # dprint('num genes greatest in curPuck', np.sum(curGreaterIdxs))

    maxExpr[curGreaterIdxs] = curPuckSums[curGreaterIdxs]
    maxExprPuck[curGreaterIdxs] = pid

    pgroup_genes = root.create_group('genes', overwrite=True)
    # pgroup_genes = zarr.empty(5, dtype=object, object_codec=numcodecs.JSON())
    # pgroup_genes_arr = root.zeros('genes', shape=(nGenes), dtype='S6')
    # pgroup_genes = root.zeros('X', shape=(1,nGenes), dtype='S6')
    pgroup_genes_arr = pgroup_genes.zeros('X', shape=(nGenes), dtype='object', object_codec=numcodecs.JSON())
    # pgroup_genes_arr[0, :] = np.asarray(gene_names)
    pgroup_genes_arr[:] = np.asarray(gene_names)

    if pid==np.max(pids):
        dprint(np.max(globalGroupX), np.max(globalGroupXout)) #  max of global nz beads in and out of regions
        dprint(np.min(globalGroupX), np.min(globalGroupXout)) #  min of global nz beads in and out of regions
        dprint(np.max(gzVal), np.max(gzValOut)) # max of global all beads in and out of all regions
        dprint('all_region_nz_array max', np.max(all_region_nz_array), tree_nodes_info[rid_to_idx_map[997]]['num_beads'])
        # normalizing
        globalGroupX[:] = np.transpose(np.nan_to_num(np.transpose(globalGroupX)/gzVal, posinf=0))[:,:] # posinf=0: if denominator (total beads) is zero, then expression must be zero too

        globalGroupXout[:] = np.transpose(np.nan_to_num(np.transpose(globalGroupXout)/gzValOut, posinf=0))[:,:] # posinf=0: if denominator (total beads) is zero, then expression must be zero too

        dprint(np.max(globalGroupX), np.max(globalGroupXout))

        # storing normalizing denominators
        gzValX = globalGroup.zeros('gzVal', shape=(nRids), chunks=(nRids), dtype='f4')
        gzValX[:] = gzVal[:]

        gzValOutX= globalGroup.zeros('gzValOut', shape=(nRids), chunks=(nRids), dtype='f4')
        gzValOutX[:] = gzValOut[:]

        # dprint(gzVal, gzValOut)



        # write out meta info
        gene_names_dict = {"data": list(gene_names), "maxExprPuck": list(maxExprPuck)}
        dprint(gene_names)
        # Serializing json
        json_object = json.dumps(gene_names_dict, indent=4)

        gene_names_file = f'{op_path}/gene_info.json'
        # Writing to sample.json
        with open(gene_names_file, "w") as outfile:
            outfile.write(json_object)





# adata = ann.AnnData(X=np.zeros((nRids, nGenes)))
# dprint(adata)


# ouptut to zarr
# zarr_filename = f'{op_path}/nz_aggr.zarr'
# adata.write_zarr(zarr_filename, (1, nGenes)) # storing 1 gene per chunk
# zfile = zarr.open(zarr_filename, mode='r+') # no need to close explicitly - https://zarr.readthedocs.io/en/stable/tutorial.html

# p1 = root.create_group('p1')
# p1X = p1.zeros('X', shape=(nRids, nGenes), chunks=(1, nGenes), dtype='i4')
# p1X[:] = 1
# p3 = root.create_group('p3', overwrite=True)
# p3X = p3.zeros('X', shape=(nRids, nGenes), chunks=(1, nGenes), dtype='f4')
# p3X[:] = 0
# p3X = zarr.open('zarr_filename/p3', mode='w', shape=(nRids, nGenes), chunks=(1, nGenes), dtype='i4')

ztest = zarr.open(zarr_filename, mode='r')
dprint(ztest.tree())
# dprint(np.array(ztest['p1/X'][1326,:]))
dprint(np.array(ztest['rids/X'])[0,3])
# dprint(ztest['genes/X'][0, :])
dprint(ztest['genes/X'][:])

exit(0)
# dprint("num leaves:", len([x for x in tree.leaves()]))

# # The file should appear in the reference space key directory
# os.listdir(reference_space_key)
#

# dprint(os.listdir(reference_space_key))

# # dprint(tree.get_structures_by_name(['Cerebrum', 'Hippocampal formation']))
# cerebrum_node=tree.get_structures_by_name(['Cerebrum'])
# hippo_form_node = tree.get_structures_by_name(['Hippocampal formation'])
# induseum_node = tree.get_structures_by_name(['Induseum griseum'])
# ammons_node = tree.get_structures_by_name(['Ammon\'s horn'])
# hippo_region_node = tree.get_structures_by_name(['Hippocampal region'])

# dprint(induseum_node)
# print("")
# print("")
# dprint(hippo_form_node)
# print("")
# print("")
# dprint(ammons_node)
# print("")
# print("")
# dprint(hippo_region_node)
# basic_groups_node=tree.get_structures_by_name(['Basic cell groups and regions'])

# nrrd_path = f'{reference_space_key}/annotation_10.nrrd'
nrrd_path = f'{reference_space_key}/annotation_25.nrrd'
readdata, header = nrrd.read(nrrd_path)
dprint(readdata.shape)
# dprint(header)
uniq  = np.unique(readdata)
dprint("uniq:", len(uniq))


# steps:
# - create tree of type from treelib import Tree  using ancestors list - use add_nodes function from analysis/s9f_dendro.py
# - traverse recursively the tree starting from root and followed by children - see get_child_info in s9f_dendro.py and note tree_nodes variable for accumulation ; initialize gene nz_counts array for leaf nodes (use R script's output for obtaining these vals)
# - on returns, add up the counts for the ancestors
# - write out json with keys as region ids mapped to array of length num genes; also regionwise bead count data available at s9_analysis/s9d/interim/aggr_num_beads_xxx.csv; global counts to compute percent outside region can be found from root's info

# what happens when multiple regions selected?? - deal in client side
