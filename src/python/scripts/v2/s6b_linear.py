"""Maps from Chuck space image coordinates to Allen atlas image coords

Usage:

python s6b_linear.py \
    nissl ids to process \
    path to input csv files with chuck space image coordinates \
    quick nii json file with transform parameters \
    path to output.csv files

Usage example:

python src/python/scripts/v2/s6b_linear.py \
    143 \
    /Users/mraj/Desktop/forgif/chuck_space_img_coords \
    /Users/mraj/Desktop/transformed_hz_png/March23_Full_CCF.json \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_allen_coords

Created by Mukund on 2022-03-29

References:

- 25 micro meter resolution
- dims = [528 320 456];
- http://dirsig.cis.rit.edu/docs/new/affine.html
- https://medium.com/geekculture/right-and-left-matrix-multiplication-d21947f195d8

"""

from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import csv
import numpy as np
import json

nissl_id = int(sys.argv[1])
ip_coords_folder = sys.argv[2]
quicknii_json_file = sys.argv[3]
op_folder = sys.argv[4]

if (nissl_id<0):
    print("mapping all")
    # nissl_ids = [141, 143]
    # nissl_ids = [123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143]
    nissl_ids = [1, 3, 133, 135, 137, 139, 143]
    nissl_ids = [i for i in np.arange(1,228,2)]
    nissl_ids.remove(5)
    nissl_ids.remove(77)
    nissl_ids.remove(167)
else:
    nissl_ids = [nissl_id]

# read json
f = open(quicknii_json_file)

# returns JSON object as
# a dictionary
qnii_data_unsorted = json.load(f)["slices"]
qnii_data = {}
img_dims = {}

for item in qnii_data_unsorted:
    nissl_id = int(item["filename"][4:7])
    anch = item["anchoring"]
    qnii_data[nissl_id] = np.array([[anch[3], anch[4], anch[5]],
                                    [anch[6], anch[7], anch[8]],
                                    [anch[0], anch[1], anch[2]]])
    img_dims[nissl_id] = {"height": item["height"], "width":item["width"]}


for nissl_id in nissl_ids:
    print (nissl_id)

    # read bead positions
    nis_idx = str(nissl_id).zfill(3)
    coords_file = ip_coords_folder+"/chuck_sp_img_coords_"+str(nis_idx)+".csv"
    input_pts = []
    with open(coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [row[0], row[1]]
            input_pts.append(point)

    pts = np.array(input_pts)

    # create normalized, homogeneous coords
    # pts = pts/np.array
    print ("width: ", img_dims[nissl_id]["width"], "height:", img_dims[nissl_id]["height"])
    den = np.array([img_dims[nissl_id]["width"], img_dims[nissl_id]["height"], 1]).reshape((3,1))
    N = len(pts)
    ones = np.ones(N).reshape((N,1))
    pts = np.concatenate((pts, ones), axis=1)
    # print(np.shape(ones),np.shape(pts))
    pts = pts.T/den
    pts = pts.T
    print("Pts normalized, homogeous:\n")
    print(pts)
    print("\n")
    print(np.amax(pts, axis=0), np.amin(pts, axis=0))

    pts = pts@qnii_data[nissl_id]
    print("\nPts in physical space:\n")
    print(pts)
    pts_physical = pts
    print("\n")
    print(np.amax(pts, axis=0), np.amin(pts, axis=0))

    # transform to allen space
    ones = np.ones(N).reshape((N,1))
    pts = np.concatenate((pts, ones), axis=1)

    print("\nPts in physical space homogenized:\n")
    print(pts)
    T_allen = np.array([[0, 0, 25, 0],
                        [-25, 0, 0, 0],
                        [0, -25, 0, 0],
                        [13175, 7975, 0, 1]])

    T_allen2 = np.array([[0.04, 0, 0, 0],
                        [0, 0.04, 0, 0],
                        [0, 0, 0.04, 0],
                        [0, 0, 0, 1]])
    pts = pts@T_allen
    pts = pts@T_allen2 # to convert to Allen image space

    print("\nPts in allen img space:\n")
    print(pts)

    print("\n")
    print(np.amax(pts, axis=0), np.amin(pts, axis=0))

    pts = list(pts)
    # write out coords after converting to int
    op_file = op_folder+"/allen_img_coords_"+nis_idx+".csv"
    with open(op_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, pt in enumerate(pts_physical):
            row = [idx, int(pt[0]), int(pt[1]), int(pt[2])]
            writer.writerow(row)



