"""
Script to generate hierarchical json data for dendrogram component in viewer.

Usage:

python s9f_dendro.py \
    inp: data root
    out: path to output folder

Usage example:

python src/python/scripts/analysis/s9f_dendro.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9e/gene_jsons_s9f \

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
op_folder = data_root+sys.argv[2]

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



tree.show()
def get_child_info(tree, child_id, name_map):

    children = tree.children(child_id)
    # dprint(children)
    data = {"label":name_map[child_id], "value":child_id}
    if (len(children)>0):
        data["children"] = []
        for cur_child_node in children:
            cur_child_id = cur_child_node.identifier
            info = get_child_info(tree, cur_child_id, name_map)
            data["children"].append(info)

    return data


data = get_child_info(tree, ancestors_list[0][0], name_map)
dprint(data)
# sub_t = tree.subtree(843)
# dprint(len(list(sub_t.expand_tree())))

op_folder = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9f"
# write out outdata as json
op_file = f'{op_folder}/regions.json'
with open(op_file, 'w') as outfile:
    json.dump(data, outfile, separators=(',', ':'))

dprint('done')
