"""
Quick script to generate mapping for ccf coords to region for Evan experiments.
Generates the coords to regionid mapping and regionid to region name mapping.

Usages:

python gen_coord_to_ccf.py \
    i/o: data root
    out: path to output folder

Example:

python src/python/scripts/temp/gen_coord_to_ccf.py \
    ~/Desktop/work/data/mouse_atlas \
    /scratch/coords_to_ccfregion \

Created by Mukund on 2022-11-28
"""


import numpy as np
import nrrd
from produtils import dprint
import sys
import csv
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import json

data_root = sys.argv[1]
op_folder = data_root+sys.argv[2]

nrrd_path = 'annotation/ccf_2017/annotation_25.nrrd'
readdata, header = nrrd.read(nrrd_path)


m,n,p = readdata.shape

## uncomment following block to generate coords to regionid mapping
# with open (op_folder+'/coords_to_ccfregion.csv', 'w') as f:
#     for x in range(m):
#         for y in range(n):
#             for z in range(p):
#                 print(x,y,z, readdata[x,y,z])
#                 f.write(str(x)+','+str(y)+','+str(z)+','+str(readdata[x,y,z])+'\n')


# for name map

reference_space_key = 'annotation/ccf_2017'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1) 

id_to_name_map = tree.get_name_map()

# write out id_to_name_map to tsv file
with open(op_folder+'/id_to_name_map.tsv', 'w') as f:
    for key in id_to_name_map.keys():
        name = id_to_name_map[key]
        structure = tree.get_structures_by_name([name])
        allen_id = structure[0]['acronym']
        f.write(str(key)+'\t'+allen_id+'\t'+name+'\n')
        # f.write(str(key)+','+allen_id+'\n')

dprint('done')
