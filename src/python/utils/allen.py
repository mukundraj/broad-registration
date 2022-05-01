import nrrd

"""
Helper functions related to allen sdk

Created by Mukund on 2022-04-27
"""

from allensdk.core.reference_space_cache import ReferenceSpaceCache

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
