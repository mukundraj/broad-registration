import nrrd

"""
Helper functions related to allen sdk

Created by Mukund on 2022-04-27
"""

from produtils import dprint
from allensdk.core.reference_space_cache import ReferenceSpaceCache

def get_cortex_layer_and_hippo_ids_lists():
    """
    Gets a list of lists with each list having region ids of regions in Allen 
    reference atlas that belong to a layer in cortex. Initially created for 
    purpose of QC heatmaps

    """

    reference_space_key = 'annotation/ccf_2017'
    resolution = 25
    rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
    # ID 1 is the adult mouse structure graph
    tree = rspc.get_structure_tree(structure_graph_id=1) 

    name_map = tree.get_name_map()
    id_map = {v: k for k, v in name_map.items()}
    hippo_id = tree.get_structures_by_name(['Hippocampal formation'])[0]['id']

    l1= [value for key, value in id_map.items() if 'layer 1' in key.lower()]
    l2l3= [value for key, value in id_map.items() if 'layer 2/3' in key.lower()]
    l4 = [value for key, value in id_map.items() if 'layer 4' in key.lower()]
    l5 = [value for key, value in id_map.items() if 'layer 5' in key.lower()]
    l6 = [value for key, value in id_map.items() if 'layer 6' in key.lower()]
    hippo = [value for key, value in id_map.items() if tree.structure_descends_from(value, hippo_id)]

    dprint(len(l1), len(l2l3), len(l4), len(l5), len(l6), len(hippo))
    # dprint(len(id_map.keys()))

    cortex_layer_and_hippo_ids = [l1, l2l3, l4,l5, l6, hippo]
    return cortex_layer_and_hippo_ids


def sample_allen_annotation(in_pts, nrrd_path):
    """
    Given a list of pts, and path to annotations nrrd, returns a list of
    sampled annotation values.

    """

    readdata, header = nrrd.read(nrrd_path)
    print("nrrd dims: ", readdata.shape)

    allen_annos = []
    # reference_space_key = 'annotation/ccf_2017/' #fixme
    # resolution = 25
    # rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
    # # ID 1 is the adult mouse structure graph
    # tree = rspc.get_structure_tree(structure_graph_id=1) 
    # name_map = tree.get_name_map()

    for pt in in_pts:
        x = int(pt[0])
        y = int(pt[1])
        z = int(pt[2])
        if (x>=readdata.shape[0] or y>= readdata.shape[1] or z>= readdata.shape[2] or \
                x<0 or y<0 or z<0):
            allen_annos.append(0) 
        else:
            id = readdata[x, y, z]
            allen_annos.append(id)

    return allen_annos


def get_allen_anno_from_id(id):

    reference_space_key = 'annotation/ccf_2017/'
    resolution = 25
    rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
    # ID 1 is the adult mouse structure graph
    tree = rspc.get_structure_tree(structure_graph_id=1) 

    code = 2
    name = 3


    # print(tree.structure_descends_from(id, cerebrum_node[0]['id']))
    # print(tree.structure_descends_from(id, basic_groups_node[0]['id']))

    return [code, name]
