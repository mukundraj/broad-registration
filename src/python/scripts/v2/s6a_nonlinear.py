"""
Use VisuAlign fiducials to convert transform image and coordinates

Usage:

python s6a_nonlinear.py \
    inp: nissl ids to process \
    inp: path to bead coords to be transformed \
    inp: path to slicer fiducials folder \
    out: path for storing intrim chuck space coords with sign change \
    out: path to output transformed bead coords

Usage example:

python src/python/scripts/v2/s6a_nonlinear.py \
    -1 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_fiducial_csvs \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_intrim_for_nlalign \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_nlaligned_coords_allbds

References:

- https://itk.org/ITKExamples/src/Filtering/ImageGrid/ResampleAnImage/Documentation.html

Created by Mukund on 2022-03-30

"""


from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.render2d as render2d

import csv
import numpy as np
import json
import subprocess
from src.python.utils.misc import get_ofixed_nis_idx

nissl_id = int(sys.argv[1])
chuck_space_coords_folder = sys.argv[2]
tfm_folder = sys.argv[3]
intrim_pts_folder = sys.argv[4]
op_nlaligned_coords_folder = sys.argv[5]

if (nissl_id<0):
    print("mapping all")
    # nissl_ids = [141, 143]
    # nissl_ids = [123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143]
    nissl_ids = [1, 3, 133, 135, 137, 139, 143]
    nissl_ids = [i for i in np.arange(1,214,2)]
    nissl_ids.remove(5)
    nissl_ids.remove(77)
    nissl_ids.remove(167)

else:
    nissl_ids = [nissl_id]

print (nissl_id)

for nissl_id in nissl_ids:
    # warped_coords
    nis_idx = str(nissl_id).zfill(3)

    intrim_pts_file = chuck_space_coords_folder+"/chuck_sp_img_coords_"+nis_idx+".csv" # pts to be transformed aka input pts

    # read input pts and convert to negative
    in_pts = []
    with open(intrim_pts_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [-row[0], -row[1]]
            in_pts.append(point)

    # write back negative
    intrim_pts_file = f'{intrim_pts_folder}/chuck_space_coords_neg_{nis_idx}.csv'
    with open(intrim_pts_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for row in in_pts:
            writer.writerow(row)

    cnissl_id = get_ofixed_nis_idx(nissl_id)
    from_fiducials_file = tfm_folder+"/"+str(cnissl_id)+"_f.csv"
    to_fiducials_file = tfm_folder+"/"+str(cnissl_id)+"_t.csv"
    # nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
    temp_output_csv_file = op_nlaligned_coords_folder+"/qnii_coords_"+nis_idx+".csv"
    # call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
    subprocess.run(["./build/cmapper2",
                    from_fiducials_file,
                    to_fiducials_file,
                    intrim_pts_file,
                    "None",
                    nis_idx,
                    temp_output_csv_file])

    warped_coords = []
    with open(temp_output_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [-row[0], -row[1]]
            warped_coords.append(point)

    with open(temp_output_csv_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for row in warped_coords:
            writer.writerow(row)
    # convert to image coords to positive
    # nl_tfmed_img_coords = [[pt[0], pt[1]] for pt in warped_coords]

    # input_pts = []
    # with open(intrim_pts_file, newline='\n') as csvfile:
    #     reader = csv.reader(csvfile)
    #     for row in reader:
    #         # print(row['first_name'], row['last_name'])
    #         row = list(map(float, row))
    #         point = [row[0], row[1]]
    #         input_pts.append(point)

    # call drawer function - pass image dims, input pts,  and output path
    # width = 4096
    # height = 3606
    # op_img_path = f'{op_nlaligned_imgs_folder}/nlaligned_{nis_idx}.png'
    # render2d.draw_image(warped_coords, width, height, op_img_path)
