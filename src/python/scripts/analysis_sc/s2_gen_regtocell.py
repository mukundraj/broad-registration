"""Create a dict mapping region id to array of celltype ids.

Usage:

python s2_gen_regtocell.py
    i/o: data root
    inp: regions tree json file
    inp: input folder with cell x bead anndata files
    inp: path to 3D nrrd file with CCF region ids
    out: path to output folder

Example:

python src/python/scripts/analysis_sc/s2_gen_regtocell.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9f/regions.json \
    /cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW \
    annotation/ccf_2017/annotation_25.nrrd \
    /single_cell/s2/s2_regtocell

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s2/s2_regtocell gs://bcdportaldata/singlecell_data/s2/s2_regtocell

Created by Mukund on 2022-11-27

"""

from produtils import dprint
import json
import sys
import anndata as ann
from scipy.sparse import csr_matrix, hstack
import numpy as np
import nrrd


data_root = sys.argv[1]
region_tree_file = data_root+sys.argv[2]
ip_folder_cxb = data_root+sys.argv[3]
nrrd_path = sys.argv[4]
readdata, header = nrrd.read(nrrd_path)
op_folder = data_root+sys.argv[5]

with open (region_tree_file, 'r') as f:
    tree = json.load(f)


nodes_list = []
def traverse(node, nodes_list):
    nodes_list.append(node)
    if 'children' in node:
        children = node['children']
        for cur_child_node in children:
            traverse(cur_child_node, nodes_list)

    return


# dprint(tree['children'])
traverse(tree, nodes_list)
dprint(len(nodes_list))


# initialize region_to_celltype with empty sets for each regionid
region_to_celltype = {}

dprint(nodes_list)
for node in nodes_list:
    region_to_celltype[node['value']] = set()

dprint(len(region_to_celltype.keys()))

# dprint(region_to_celltype)
# dprint('num regions', len(region_to_celltype.keys()), ctr)

# iterate over pucks
# get puck ids
start_pid = 1
end_pid = 207
pids = list(range(int(start_pid), int(end_pid)+1, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

cells = None
for pids_idx, pid in enumerate(pids):
    assert(pid!=5 and pid!=77 and pid!=167)
    nis_id_str = str(pid).zfill(3)
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # iterate over cellidx
    ip_counts_file  = f'{ip_folder_cxb}/dd{str(apid)}_CTMapping.h5ad'
    dprint('pid', pid)
    counts = ann.read_h5ad(ip_counts_file)
    positions = counts.var[['CCF3D_x', 'CCF3D_y', 'CCF3D_z']].values
    position_inds = counts.var.rowName.values.astype(int)
    counts_X = csr_matrix(counts.X).transpose()
    dprint(position_inds.shape, counts_X.shape, positions.shape)


    cells = list(counts.obs.index)
    cells = [x.split('=')[1] for x in cells]
    for cell_idx, cell in enumerate(cells):
        # get cell_scores for for all beads for this cell
        specific_cell_scores = counts_X.getcol(cell_idx)
        spec_cell_scores_dense = np.squeeze(np.array(specific_cell_scores.todense())).astype(float)
        # dprint(spec_cell_scores_dense.shape)


        # get inds of beads > 0.3
        bead_inds = np.where(spec_cell_scores_dense > 0.3)[0]
        # dprint(bead_inds.shape)
        # dprint(bead_inds)


        # for chosen beads get region ids and region_to_celltype[region_id].add(cellidx)
        for bead_idx in bead_inds:
            bead_pos = positions[bead_idx,:]
            x = int(bead_pos[0])
            y = int(bead_pos[1])
            z = int(bead_pos[2])
            if (x>=readdata.shape[0] or y>= readdata.shape[1] or z>= readdata.shape[2] or \
                    x<0 or y<0 or z<0):
                id = 0
            else:
                id = readdata[x, y, z]
            region_to_celltype[id].add(cell_idx)




# hydrate sets of nonroot elements with their children's sets

def hydrate(node):
    if 'children' in node:
        children = node['children']
        for cur_child_node in children:
            hydrate(cur_child_node)
            if 'celltype_ids' in node:
                node['celltype_ids'] = node['celltype_ids'].union(cur_child_node['celltype_ids'])
            else:
                node['celltype_ids'] = cur_child_node['celltype_ids']
    else:
        node['celltype_ids'] = region_to_celltype[node['value']]
    return


hydrate(tree)

nodes_list = []
traverse(tree, nodes_list)

for node in nodes_list:
    region_to_celltype[node['value']] = node['celltype_ids']

# convert sets to lists
for key in region_to_celltype.keys():
    region_to_celltype[key] = list(region_to_celltype[key])

# write out region_to_celltype as json with no space
with open(f'{op_folder}/region_to_celltype.json', 'w',) as f:
    json.dump(region_to_celltype, f, separators=(',', ':'))

cells_to_idx = {}
for cell_idx, cell in enumerate(cells):
    cells_to_idx[cell] = cell_idx

# write out cells_to_idx as json with no space
with open(f'{op_folder}/mappedCellType_to_idx.json', 'w',) as f: # mappedCellTypes as opposed to all cellTypes identified during clustering
    json.dump(cells_to_idx, f, separators=(',', ':'))

