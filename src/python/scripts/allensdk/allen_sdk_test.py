"""
Script to test Allen SDK

Created by Mukund on 2022-04-26


Usage example:

python src/python/scripts/allensdk/allen_sdk_test.py

References:
- https://allensdk.readthedocs.io/en/latest/reference_space.html
- https://allensdk.readthedocs.io/en/latest/_static/examples/nb/reference_space.html#Downloading-an-annotation-volume
"""


import os
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import nrrd

reference_space_key = 'annotation/ccf_2017/'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1)


# annotation, meta = rspc.get_annotation_volume()

# # The file should appear in the reference space key directory
# os.listdir(reference_space_key)
#

print(os.listdir(reference_space_key))

print(tree.get_structures_by_name(['Cerebrum']))
cerebrum_node=tree.get_structures_by_name(['Cerebrum'])
basic_groups_node=tree.get_structures_by_name(['Basic cell groups and regions'])

nrrd_path = f'{reference_space_key}/annotation_25.nrrd'
readdata, header = nrrd.read(nrrd_path)
print(readdata.shape)
print(header)

print (readdata[268, 160, 228])

id = readdata[268, 160, 228]

name_map = tree.get_name_map()
print(id, name_map[id])
print(id, cerebrum_node[0]['id'])

print(tree.structure_descends_from(id, cerebrum_node[0]['id']))
print(tree.structure_descends_from(id, basic_groups_node[0]['id']))
# labels = [readdata[pt[0]][pt[1]] for pt in nissl_img_corrds]
#


