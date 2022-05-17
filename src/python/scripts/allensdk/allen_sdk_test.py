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
from produtils import dprint
import numpy as np

reference_space_key = 'annotation/ccf_2017/'
resolution = 25
# resolution = 10
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# rspc = ReferenceSpaceCache(resolution, reference_space_key)
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1)

annotation, meta = rspc.get_annotation_volume()

# # The file should appear in the reference space key directory
# os.listdir(reference_space_key)
#

dprint(os.listdir(reference_space_key))

# dprint(tree.get_structures_by_name(['Cerebrum', 'Hippocampal formation']))
cerebrum_node=tree.get_structures_by_name(['Cerebrum'])
hippo_form_node = tree.get_structures_by_name(['Hippocampal formation'])
induseum_node = tree.get_structures_by_name(['Induseum griseum'])
ammons_node = tree.get_structures_by_name(['Ammon\'s horn'])
hippo_region_node = tree.get_structures_by_name(['Hippocampal region'])

dprint(induseum_node)
print("")
print("")
dprint(hippo_form_node)
print("")
print("")
dprint(ammons_node)
print("")
print("")
dprint(hippo_region_node)
basic_groups_node=tree.get_structures_by_name(['Basic cell groups and regions'])

nrrd_path = f'{reference_space_key}/annotation_25.nrrd'
readdata, header = nrrd.read(nrrd_path)
dprint(readdata.shape)
dprint(header)

unique, counts = np.unique(readdata, return_counts=True)
count_dict = dict(zip(unique, counts))
dprint(len(count_dict.keys()))
dprint(count_dict[1089])
dprint(count_dict[19])

dprint(hippo_region_node[0]['id'])
tofind = hippo_region_node[0]['id']
if (tofind in readdata):
    dprint("found it", tofind)
else:
    dprint("not found")


dprint (readdata[268, 160, 228])

id = readdata[268, 160, 228]

name_map = tree.get_name_map()
dprint(id, name_map[id])
dprint(id, cerebrum_node[0]['id'])

dprint(tree.structure_descends_from(id, cerebrum_node[0]['id']))
dprint(tree.structure_descends_from(id, basic_groups_node[0]['id']))
# labels = [readdata[pt[0]][pt[1]] for pt in nissl_img_corrds]
#


