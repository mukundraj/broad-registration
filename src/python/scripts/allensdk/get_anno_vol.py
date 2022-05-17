"""
Script to get Allen SDK annotation volume

Created by Mukund on 2022-02-28

Usage example:

python src/python/scripts/allensdk/get_anno_vol.py

References:
- https://allensdk.readthedocs.io/en/latest/reference_space.html
- https://allensdk.readthedocs.io/en/latest/_static/examples/nb/reference_space.html#Downloading-an-annotation-volume
"""


import os
from allensdk.core.reference_space_cache import ReferenceSpaceCache

reference_space_key = 'annotation/ccf_2017'
resolution = 25
resolution = 10
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1) 

annotation, meta = rspc.get_annotation_volume()

# The file should appear in the reference space key directory
os.listdir(reference_space_key)
