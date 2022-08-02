"""
Script to generate hierarchical json data for dendrogram component in viewer.

Usage:

python s9f_dendro.py \

Usage example:

python src/python/scripts/analysis/s9f_dendro.py \

Created by Mukund on 2022-08-02

"""

import os
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import nrrd
from produtils import dprint
import numpy as np
from treelib import Node, Tree 

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

# tree.show()
# sub_t = tree.subtree(843)
# dprint(len(list(sub_t.expand_tree())))

dprint('done')
