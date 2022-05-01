"""
Reads in allen bead positions and outputs table with CCF regions

Usage:

python bead_to_ccf_labels.py \
    inp: nissl_id
    inp: path to annotation nrrd file
    inp: path to folder with allen coords
    out: path to folder with allen region info table
    out: path to folder with allen region rgb info only (for smaller file size)

Usage example:

python src/python/scripts/allensdk/bead_to_ccf_labels.py \
    -1 \
    annotation/ccf_2017/annotation_25.nrrd \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_allen_coords \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/bead_ccf_rgb

References:

- 25 micro meter resolution
- dims = [528 320 456] in PIR coordinate system (as opposed to RAS+ of NIfTI)

Created by Mukund on 2022-04-27
"""

from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.allen as allen
import csv
import numpy as np
import math
import subprocess
import nrrd
from allensdk.core.reference_space_cache import ReferenceSpaceCache

nissl_id = int(sys.argv[1])
nrrd_path = sys.argv[2]
allen_coords_folder = sys.argv[3]
op_folder_anno = sys.argv[4]
op_folder_rgb = sys.argv[5]

# read globals
if (nissl_id<0):
    print("mapping all")
    nissl_ids = [i for i in np.arange(1,208,2)]
    nissl_ids.remove(5)
    nissl_ids.remove(77)
    nissl_ids.remove(167)
    # nissl_ids = [97, 143]
else:
    nissl_ids = [nissl_id]

readdata, header = nrrd.read(nrrd_path)
print("nrrd dims: ", readdata.shape)
# print(header)
reference_space_key = 'annotation/ccf_2017/'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1)
name_map = tree.get_name_map()

for nissl_id in nissl_ids:
    print(f'processing {nissl_id}')
    nis_id_str = str(nissl_id).zfill(3)
    allen_coods_file = f'{allen_coords_folder}/allen_img_coords_{nis_id_str}.csv'
    in_beads = []
    in_beads_names = []
    in_beads_ids = []
    with open(allen_coods_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            x = int(row[1])
            y = int(row[2])
            z = int(row[3])
            if (x>=readdata.shape[0] or y>= readdata.shape[1] or z>= readdata.shape[2] or \
                    x<0 or y<0 or z<0):
                continue
            id = readdata[x, y, z]
            if (id>0):
                name = name_map[id]
                # structure = tree.get_structures_by_name([name])
                # allen_id = structure[0]['acronym']
                # print(row['first_name'], row['last_name'])
                # [code, name] = allen.get_allen_anno_from_id(id)
                # row = list(map(float, row))
                # point = [-row[0], -row[1]]
                # in_pts.append(point)
                in_beads.append([x,y,z])
                in_beads_names.append(name)
                in_beads_ids.append(id)

    in_beads_structures = tree.get_structures_by_name(in_beads_names)
    in_beads_acronyms = list(map(lambda x: x["acronym"], in_beads_structures))
    in_beads_rgb = list(map(lambda x: x["rgb_triplet"], in_beads_structures))
    in_beads_babylon = []
    # print(in_beads_acronyms)
    # print(len(in_beads))

    op_file = f'{op_folder_anno}/allen_anno_data_{nis_id_str}.csv'
    with open(op_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, row in enumerate(in_beads):
            line = [row[0], row[1], row[2], in_beads_acronyms[idx], in_beads_names[idx], in_beads_rgb[idx][0], in_beads_rgb[idx][1], in_beads_rgb[idx][2]]
            writer.writerow(line)
            in_beads_babylon.append([row[0], 320-row[1], 456-row[2]])

    # # To invert the y axis for visualization in Babylon js coordinates
    # T_babylon = np.array([[1, 0, 0, 0],
    #                     [0, -1, 0, 0],
    #                     [0, 0, 1, 0],
    #                     [0, 320, 0, 1]])
    # in_beads = np.array(in_beads)

    # in_beads_babylon = in_beads.T@T_babylon

    op_file = f'{op_folder_rgb}/babylon_coords_rgb_{nis_id_str}.csv'
    with open(op_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, row in enumerate(in_beads_babylon):
            line = [row[0], row[1], row[2], round(float(in_beads_rgb[idx][0])/255,2), round(float(in_beads_rgb[idx][1])/255, 2) , round(float(in_beads_rgb[idx][2])/255,2)]
            writer.writerow(line)

