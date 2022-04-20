"""
For coordinate normalization after accounting for ggplot induced padding

Usage:

python s4a_coord_norm.py \
    inp: nissl_id \
    inp: path to corners csv file \
    inp: path to folder with bead coord files \
    out: path to folder to place normalized coord files

Usage example:

python src/python/scripts/v2/s4a_coord_norm.py \
    -1 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/corners.csv \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_raw_data/bead_coords \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s2_seg_ids/normalized_coords


Created by Mukund on 2022-03-18

References:

https://stackoverflow.com/questions/57991865/numpy-finding-the-tranformation-matrix-given-2-sets-of-4-points-in-3d-euclidia

"""

from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import csv
import numpy as np
import src.python.utils.io as io
import sys

nissl_id = int(sys.argv[1])
corners_file = sys.argv[2]
bead_coords_folder = sys.argv[3]
op_folder = sys.argv[4]

# read globals
if (nissl_id<0):
    print("mapping all")
    # nissl_ids = [141, 143]
    # nissl_ids = [123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143]
    nissl_ids = [1, 3, 133, 135, 137, 139, 143]
    nissl_ids = [i for i in np.arange(1,228,2)]
    nissl_ids.remove(5)
    nissl_ids.remove(77)
    nissl_ids.remove(167)
    nissl_ids = [65, 67]
else:
    nissl_ids = [nissl_id]

# read corners csv
corners = io.get_corners_info(corners_file)

for nissl_id in nissl_ids:
    print("Mapping nissl_id:", nissl_id)


    topleft_x = corners[nissl_id]['topleft_x']
    topleft_y = corners[nissl_id]['topleft_y']

    botright_x = corners[nissl_id]['botright_x']
    botright_y = corners[nissl_id]['botright_y']

    topright_x = corners[nissl_id]['topright_x']
    topright_y = corners[nissl_id]['topright_y']

    botleft_x = corners[nissl_id]['botleft_x']
    botleft_y = corners[nissl_id]['botleft_y']

    nis_id_str_nozfill = str(nissl_id)
    nis_id_str = str(nissl_id).zfill(3)
    ip_file = f'{bead_coords_folder}/coords_{nis_id_str_nozfill}.csv'
    op_file = f'{op_folder}/normed_coords_{nis_id_str}.csv'

    print(topleft_x, topleft_y, botright_x, botright_y, topright_x, topright_y, botleft_x, botleft_y)
    print(ip_file, op_file )


    s = np.array([[topleft_x, topleft_y, 1], 
                  [topright_x, topright_y, 1],
                  [botleft_x, botleft_y, 1],
                  [botright_x, botright_y, 1]], dtype=float).T

    d = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float).T

    M_rec,resid,rank,sing = np.linalg.lstsq(s.T,d.T)
    M_rec = M_rec.T


    t = np.array([[671, 575, 1], [2566, 2796, 1]],dtype=float)
    t = t.T

    # read input
    input_pts = []
    with open(ip_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [row[0], row[1]]
            input_pts.append(point)

    print("Number of positions read: ", len(input_pts))

    # format input
    N = len(input_pts)
    pts = np.array(input_pts)
    arrays = [pts, np.ones(N).reshape(N,1)]
    pts = np.hstack((pts, np.ones(N).reshape(N,1))).T

    # transform input
    tfmed_pts = M_rec@pts

    tfmed_pts = list(tfmed_pts.T)
    # print(tfmed_pts)

    # write output
    with open(op_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for row in tfmed_pts:
            writer.writerow([round(row[0], 5), round(row[1], 5)])
