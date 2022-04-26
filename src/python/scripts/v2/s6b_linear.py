"""Maps from Chuck space image coordinates to Allen atlas image coords

Usage:

python s6b_linear.py \
    inp: nissl ids to process \
    inp: path to input csv files with chuck space image coordinates \
    inp: quick nii json file with transform parameters \
    out: path to output.csv files in allen coordinates
    out: path to output.csv files visualization in babylon

Usage example(s):

python src/python/scripts/v2/s6b_linear.py \
    -1 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_nlaligned_coords  \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png_copy/QuickNii_4-17-22-mod.json \
    4096 3606 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_allen_coords \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_babylon_coords

Created by Mukund on 2022-03-29

References:

- 25 micro meter resolution
- dims = [528 320 456] in PIR coordinate system (as opposed to RAS+ of NIfTI)
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
img_width = int(sys.argv[4])
img_height = int(sys.argv[5])
op_folder_allen = sys.argv[6]
op_folder_babylon = sys.argv[7]

if (nissl_id<0):
    print("mapping all")
    # nissl_ids = [141, 143]
    # nissl_ids = [123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143]
    nissl_ids = [1, 3, 133, 135, 137, 139, 143]
    nissl_ids = [i for i in np.arange(1,228,2)]
    nissl_ids = [i for i in np.arange(1,208,2)]
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
    print(nissl_id)
    if nissl_id > 208:
        break
    qnii_data[nissl_id] = np.array([[anch[3], anch[4], anch[5]],
                                    [anch[6], anch[7], anch[8]],
                                    [anch[0], anch[1], anch[2]]])
    img_dims[nissl_id] = {"height": item["height"], "width":item["width"]}


for nissl_id in nissl_ids:
    print (nissl_id)

    # read bead positions
    nis_idx = str(nissl_id).zfill(3)
    # coords_file = ip_coords_folder+"/chuck_sp_img_coords_"+str(nis_idx)+".csv"
    coords_file = ip_coords_folder+"/qnii_coords_"+str(nis_idx)+".csv"
    input_pts = []
    with open(coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [row[0], row[1]]
            input_pts.append(point)

    # convert from left/right coordinate from left of image to left of animal
    input_pts = [[img_width-pt[0], pt[1]] for pt in input_pts]

    pts = np.array(input_pts)

    # create normalized, homogeneous coords
    # pts = pts/np.array
    print ("width: ", img_dims[nissl_id]["width"], "height:", img_dims[nissl_id]["height"])
    den = np.array([img_dims[nissl_id]["width"], img_dims[nissl_id]["height"], 1]).reshape((3,1))
    N = len(pts)
    ones = np.ones(N).reshape((N,1))
    pts = np.concatenate((pts, ones), axis=1)
    # print(np.shape(ones),np.shape(pts))
    print("den", den)
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

    # to convert to Allen image space
    T_allen2 = np.array([[0.04, 0, 0, 0],
                        [0, 0.04, 0, 0],
                        [0, 0, 0.04, 0],
                        [0, 0, 0, 1]])

    # To invert the y axis for visualization in Babylon js coordinates
    T_allen3 = np.array([[1, 0, 0, 0],
                        [0, -1, 0, 0],
                        [0, 0, 1, 0],
                        [0, 320, 0, 1]])

    pts_allen = pts@T_allen@T_allen2

    pts_babylon = pts_allen@T_allen3

    print("\nPts in allen img space:\n")
    print(pts_allen)

    print("\n")
    print(np.amax(pts_allen, axis=0), np.amin(pts_allen, axis=0))

    pts = list(pts)

    # write out coords after converting to int
    op_file_allen = op_folder_allen+"/allen_img_coords_"+nis_idx+".csv"
    with open(op_file_allen, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, pt in enumerate(pts_allen):
            row = [idx, int(pt[0]), int(pt[1]), int(pt[2])]
            writer.writerow(row)

    op_file_babylon = op_folder_babylon+"/babylon_img_coords_"+nis_idx+".csv"
    with open(op_file_babylon, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, pt in enumerate(pts_babylon):
            row = [idx, int(pt[0]), int(pt[1]), int(pt[2])]
            writer.writerow(row)



