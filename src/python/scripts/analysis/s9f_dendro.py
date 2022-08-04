"""
Script to generate hierarchical json data for dendrogram component in viewer.

Usage:

python s9f_dendro.py \
    inp: data root
    inp: path to regionwise aggregated bead count data
    inp inp: start_pid end_pid
    out: path to output folder

Usage example:

python src/python/scripts/analysis/s9f_dendro.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/interim \
    1 207 \
    /data_v3_nissl_post_qc/s9_analysis/s9f \

Created by Mukund on 2022-08-02

"""

import os
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import nrrd
from produtils import dprint
import numpy as np
from treelib import Node, Tree 
import json
import sys

data_root = sys.argv[1]
ip_num_beads_data_folder = data_root+sys.argv[2]
start_pid = sys.argv[3]
end_pid = sys.argv[4]
op_folder = data_root+sys.argv[5]

reference_space_key = 'annotation/ccf_2017/'
resolution = 25
# resolution = 10
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# rspc = ReferenceSpaceCache(resolution, reference_space_key)
# ID 1 is the adult mouse structure graph
allentree = rspc.get_structure_tree(structure_graph_id=1)
name_map = allentree.get_name_map()

# annotation, meta = rspc.get_annotation_volume()

# dprint(len(tree.get_ancestor_id_map()))
ancestors_list = list(allentree.get_ancestor_id_map().values())


# tree.create_node("Harry", 7)  # root node
# tree.show()
# node = tree.get_node(7)
# dprint(node)

## start: creating tree in treelib format
def add_nodes(tree, nodes_list, name_map):
    nodes_list.reverse()

    if len(nodes_list)==1:
        # add root node
        tree.create_node(name_map[nodes_list[0]], nodes_list[0])

    else:

        # add non root node and children
        for idx in range(len(nodes_list)-1):

            # node = tree.get_node("harry")
            parent_id = nodes_list[idx]
            child_id = nodes_list[idx+1]

            parent_node = tree.get_node(parent_id)
            child_node = tree.get_node(child_id)

            # if parent exists add child
            if parent_node is not None and child_node == None:

                tree.create_node(name_map[child_id], child_id, parent=parent_id)

            elif child_node is not None:
                continue
            else:
                dprint("Parent doesn't exist", parent_id)
                exit(0)


tree = Tree()

dprint(len(ancestors_list))
for nodes_list in ancestors_list:
    add_nodes(tree, nodes_list, name_map)

dprint('tree depth:', tree.depth())

## end: creating tree in treelib format

## start: get puck with max beads for each region (leaf and nonleaf)

rcounts = {}

# get puck ids
pids = list(range(int(start_pid), int(end_pid)+1, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

# for each puck
for pids_idx, pid in enumerate(pids):
    assert(pid!=5 and pid!=77 and pid!=167)
    nis_id_str = str(pid).zfill(3)

    rcounts[pid] = {}
    for rid in name_map.keys():
        rcounts[pid][rid] = {"sum":0}

    # read puck data and get region ids
    aggr_counts_file = f'{ip_num_beads_data_folder}/aggr_num_beads_{nis_id_str}.csv'
    bead_counts = np.genfromtxt(aggr_counts_file, delimiter=',', skip_header=1)
    m,n = np.shape(bead_counts)
    dprint(m, n)

    # for each region in puck
    for puck_reg_idx in range(m):
        puck_reg = int(bead_counts[puck_reg_idx, 0])
        puck_reg_count = int(bead_counts[puck_reg_idx,1])
        # dprint(puck_reg, pid, puck_reg_idx, puck_reg_count)
        # for each possible region
        dprint(puck_reg, pid, puck_reg_idx, puck_reg_count)
        if (puck_reg!=0):
            for rid in name_map.keys():
                # check for puck region equals same as or descends from possible region
                if puck_reg_idx==rid or allentree.structure_descends_from(puck_reg, rid):

                    # if so, update values of possible region in rcounts
                    rcounts[pid][rid]["sum"] += puck_reg_count
                    # dprint(pid, rid, rcounts[pid][rid])
                    # pass

# now, for each region, identify pidx with max bead counts
rcounts_maxvals = {}
for rid in name_map.keys():
    rcounts_maxvals[rid]={"maxval":0, "maxval_pidx":-1, "maxval_pid":-1}
    for pids_idx, pid in enumerate(pids):
        assert(pid!=5 and pid!=77 and pid!=167)
        nis_id_str = str(pid).zfill(3)
        if rcounts[pid][rid]["sum"] > rcounts_maxvals[rid]["maxval"]:
            rcounts_maxvals[rid]["maxval"] = rcounts[pid][rid]["sum"]
            rcounts_maxvals[rid]["maxval_pidx"] = pids_idx
            rcounts_maxvals[rid]["maxval_pid"] = pid



# dprint(rcounts_maxvals)
## end: get puck with max beads for each region (leaf and nonleaf)

# exit(0)

## start: creation of viewer dendrogram data

tree.show()
def get_child_info(tree, child_id, name_map):

    children = tree.children(child_id)
    # dprint(children)
    if rcounts_maxvals[child_id]["maxval_pidx"] > -1:
        actions = [{"className":"action fa fa-level-up", "title":f'Jump to srno {rcounts_maxvals[child_id]["maxval_pidx"]+1} containing: {name_map[child_id]}', "maxval_pidx":rcounts_maxvals[child_id]["maxval_pidx"], "maxval_pid": rcounts_maxvals[child_id]["maxval_pid"]}]
    else:
        actions = []
    data = {"label":name_map[child_id], "value":child_id, "actions": actions }
    if (len(children)>0):
        data["children"] = []
        for cur_child_node in children:
            cur_child_id = cur_child_node.identifier
            info = get_child_info(tree, cur_child_id, name_map)
            data["children"].append(info)
            data["actions"] = actions

    return data


data = get_child_info(tree, ancestors_list[0][0], name_map)
dprint(data)

## end: creation of viewer dendrogram data

# sub_t = tree.subtree(843)
# dprint(len(list(sub_t.expand_tree())))

# start: writing out json
op_file = f'{op_folder}/regions.json'
with open(op_file, 'w') as outfile:
    json.dump(data, outfile, separators=(',', ':'))

# end: writing out json

dprint('done')
