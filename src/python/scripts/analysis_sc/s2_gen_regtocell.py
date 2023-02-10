"""Create a dict mapping region id to array of celltype ids. This is needed for
dendrogram componet in SingleCell tab.

Usage:

python s2_gen_regtocell.py
    i/o: data root
    inp: regions tree json file
    inp: input folder with cell x bead anndata files
    inp: path to 3D nrrd file with CCF region ids [REMOVE/OBSOLETE]
    inp: path toward v3 CCF labels 
    inp: path to allen name to acronym map
    out: path to output folder

Example:

python src/python/scripts/analysis_sc/s2_gen_regtocell.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9f/regions.json \
    /cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW \
    annotation/ccf_2017/annotation_25.nrrd \
    /cell_spatial/s1/cellspatial_data/cellscores_cshl \
    /v3/s1/allen_name_to_acro_map.csv \
    /single_cell/s2/s2_regtocell \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s2/s2_regtocell gs://bcdportaldata/batch_230131/singlecell_data/s2/s2_regtocell_230208

References:

https://stackoverflow.com/questions/10461531/merge-and-sum-of-two-dictionaries

Created by Mukund on 2022-11-27

"""

from produtils import dprint
import json
import sys
import anndata as ann
from scipy.sparse import csr_matrix, hstack
import numpy as np
import nrrd
import csv


data_root = sys.argv[1]
region_tree_file = data_root+sys.argv[2]
ip_folder_cxb = data_root+sys.argv[3]
nrrd_path = sys.argv[4]
readdata, header = nrrd.read(nrrd_path)
cellspatial_puck_folder = data_root+sys.argv[5]
allen_name_id_map_file = data_root+sys.argv[6]
op_folder = data_root+sys.argv[7]

with open (region_tree_file, 'r') as f:
    tree = json.load(f)


nodes_list = []
nodes_dict = {}
def traverse(node, nodes_list, nodes_dict):
    nodes_list.append(node)
    nodes_dict[node['value']] = node['label']
    if 'children' in node:
        children = node['children']
        for cur_child_node in children:
            traverse(cur_child_node, nodes_list, nodes_dict)

    return


# dprint(tree['children'])
traverse(tree, nodes_list, nodes_dict)
dprint(len(nodes_list))


# initialize region_to_celltype with empty sets for each regionid
region_to_celltype = {}

# dprint(nodes_list)
for node in nodes_list:
    # region_to_celltype[node['value']] = set()
    region_to_celltype[node['value']] = {}

# dprint(len(region_to_celltype.keys()))

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

cell_idx_totals = {} # stores the total number of cells of each celltype    

cells = None
for pids_idx, pid in enumerate(pids):
    assert(pid!=5 and pid!=77 and pid!=167)
    nis_id_str = str(pid).zfill(3)
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment


    puck_folder = f'{cellspatial_puck_folder}/puck{pid}'
    coords_csv_name = f'{puck_folder}/coords.csv'
    # coords_csv_name = f'{cellspatial_puck_folder}/coords_{pid}.csv'

    nis_id_str = str(apid).zfill(3)
    # dprint(labels_csv_file)

    # read name_to_acro_map_file csv

    region_name_to_id_map = {}
    with open(allen_name_id_map_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            region_name_to_id_map[row[2]] = int(row[0])
        region_name_to_id_map['OUT'] = 0
        region_name_to_id_map['NA'] = 0


    region_names = []
    region_ids = []
    with open(coords_csv_name, newline='\n') as csvfile:
        reader = csv.reader(csvfile, delimiter=':')
        next(reader)
        for row in reader:
            # split_row = row[0].split(':')
            region_names.append(row[2])
            region_ids.append(region_name_to_id_map[row[2]])


    dprint('region_names length', len(region_names))
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
    search_cell_idx = 989
    # dprint(f'check_cell {search_cell_idx}', cells[search_cell_idx])
    for cell_idx, cell in enumerate(cells):
        # get cell_scores for for all beads for this cell
        specific_cell_scores = counts_X.getcol(cell_idx)
        spec_cell_scores_dense = np.squeeze(np.array(specific_cell_scores.todense())).astype(float)
        # dprint(spec_cell_scores_dense.shape)


        # get inds of beads > 0.3
        bead_inds = np.where(spec_cell_scores_dense > 0.3)[0]
        # dprint(bead_inds.shape)
        # dprint(bead_inds)

        # if cell_idx==search_cell_idx:
        #     dprint('bead_inds', bead_inds)

        # for chosen beads get region ids and region_to_celltype[region_id].add(cellidx)
        for bead_idx in bead_inds:
            id = 0
            if region_names[bead_idx]!='OUT':
                # get region
                id = region_ids[bead_idx]
            if id==0:
                continue

            # if(cell_idx==search_cell_idx and id == 581): # 129 - V3, 581 - TRS
            #     dprint(f'bead_idx in RID {id}', bead_idx, spec_cell_scores_dense[bead_idx])
            # if (cell not in region_to_celltype[id]): ## check error - cell_idx instead of cell?
            if (cell_idx not in region_to_celltype[id]): ## check error - cell_idx instead of cell?
                region_to_celltype[id][cell_idx] = 1
            else:
                region_to_celltype[id][cell_idx] += 1

            if (cell_idx not in cell_idx_totals):
                cell_idx_totals[cell_idx] = 1
            else:
                cell_idx_totals[cell_idx] += 1
            # dprint(bead_idx, x, y, z, id)


# normalize celltype counts by total number of cells of that celltype
for region_id in region_to_celltype:
    for cell_idx in region_to_celltype[region_id]:
        region_to_celltype[region_id][cell_idx] /= cell_idx_totals[cell_idx]
        region_to_celltype[region_id][cell_idx] = round(region_to_celltype[region_id][cell_idx], 3)


# hydrate sets of nonroot elements with their children's sets

def hydrate(node):
    if 'children' in node:
        children = node['children']
        for cur_child_node in children:
            hydrate(cur_child_node)
            if 'celltype_ids' in node:
                # node['celltype_ids'] = node['celltype_ids'].union(cur_child_node['celltype_ids'])
                x = node['celltype_ids']
                y = cur_child_node['celltype_ids']
                node['celltype_ids'] = {k: round(x.get(k, 0) + y.get(k, 0), 3) for k in set(x) | set(y)} # merge dicts
            else:
                # node['celltype_ids'] = cur_child_node['celltype_ids'].union(region_to_celltype[node['value']])
                x = cur_child_node['celltype_ids']
                y = region_to_celltype[node['value']]
                node['celltype_ids'] = {k: round(x.get(k, 0) + y.get(k, 0), 3) for k in set(x) | set(y)} # merge dicts
    else:
        if 'celltype_ids' not in node:
            node['celltype_ids'] = region_to_celltype[node['value']]
        else:
            # node['celltype_ids'] = node['celltype_ids'].union(region_to_celltype[node['value']])
            x = node['celltype_ids']
            y = region_to_celltype[node['value']]
            node['celltype_ids'] = {k: round(x.get(k, 0) + y.get(k, 0), 3) for k in set(x) | set(y)} # merge dicts

    return


hydrate(tree)

nodes_list = []
nodes_dict = {}
traverse(tree, nodes_list, nodes_dict)

for node in nodes_list:
    region_to_celltype[node['value']] = node['celltype_ids']

# convert sets to lists
for key in region_to_celltype.keys():
    # region_to_celltype[key] = list(region_to_celltype[key])
    region_to_celltype[key] = region_to_celltype[key]

# write out region_to_celltype as json with no space
with open(f'{op_folder}/region_to_celltype.json', 'w',) as f:
    json.dump(region_to_celltype, f, separators=(',', ':'))

cells_to_idx = {}
for cell_idx, cell in enumerate(cells):
    cells_to_idx[cell] = cell_idx

# write out cells_to_idx as json with no space
with open(f'{op_folder}/mappedCellType_to_idx.json', 'w',) as f: # mappedCellTypes as opposed to all cellTypes identified during clustering
    json.dump(cells_to_idx, f, separators=(',', ':'))

