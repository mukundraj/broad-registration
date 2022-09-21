"""
Helper functions related to general data structures

Created by Mukund on 2022-09-13
"""

from treelib import Node, Tree
from produtils import dprint

### build tree from Allen ancestor_list array
def add_nodes(tree, nodes_list, name_map):
    # nodes_list.reverse()
    nodes_list = list(reversed(nodes_list))

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
