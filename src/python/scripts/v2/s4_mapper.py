""" Mapper function that reads in bead coordinates and transforms to chuck
(aligned nissl) space.
 
Usage:

python script.py\
    nis_idx (nissl slide index) or negative value to process all images\
    path to corners file \
    path to folder with bead coordinates files \
    img_width_ss_rescaled img_height_ss_rescaled \
    path to histolozee ss project xml with slide-seq initial transforms\
    img_width_ss_tfmed img_height_ss_tfmed \
    path to folder with transform files ({nis_idx}_f.csv {nis_idx}_t.csv)\
    img_width_nis_tfmed img_height_nis_tfmed \
    path to folder to store outputs \

Usage example:

python src/python/scripts/v2/s4_mapper.py \
    97 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/corners.csv \
    /Users/mraj/Desktop/work/projects/active/broad-registration/output/ss/bead_coords \
    1728 1728 \
    "/Users/mraj/Desktop/sample-hz/hz-project-ss.zee" \
    2048 2048 \
    "/Users/mraj/Desktop/forgif/transforms" \
    4096 3606 \
    /Users/mraj/Desktop/forgif/chuck_space_img_coords \

Created by Mukund on 2022-03-21

"""
from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import csv
import numpy as np
import src.python.utils.io as io
import src.python.utils.transforms as tfms
import math
import subprocess


nissl_id = int(sys.argv[1])
corners_file = sys.argv[2]
bead_coords_folder = sys.argv[3]
img_width_ss_rescaled = int(sys.argv[4])
img_height_ss_rescaled = int(sys.argv[5])
hz_project_ss_file = sys.argv[6]
img_width_ss_tfmed = int(sys.argv[7])
img_height_ss_tfmed = int(sys.argv[8])
tfm_folder = sys.argv[9]
img_width_nis_tfmed = int(sys.argv[10])
img_height_nis_tfmed = int(sys.argv[11])
op_folder = sys.argv[12]

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
    # nissl_ids = [97, 143]
else:
    nissl_ids = [nissl_id]

# read corners csv
corners = io.get_corners_info(corners_file)

for nissl_id in nissl_ids:
    print("Mapping nissl_id:", nissl_id)

    # read bead positions
    bead_coords_file = bead_coords_folder+"/coords_"+str(nissl_id)+".csv"
    input_pts = []
    with open(bead_coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [row[0], row[1]]
            input_pts.append(point)

    pts = np.array(input_pts)
    print(np.amax(pts, axis=0), np.amin(pts, axis=0))
    max_x, max_y = np.amax(pts, axis=0)
    min_x, min_y = np.amin(pts, axis=0)
    extents = {"max_x":max_x, "max_y":max_y, "min_x":min_x, "min_y":min_y}

    # input_pts.append([5409, 4798])
    x = 100
    # input_pts = [[268+10, 274+10], [5497-10, 5497-10]]
    # input_pts = []
    # input_pts.append([268+x, 274+x])
    # input_pts.append([5497-x, 5497-x])
    # input_pts.append([2882, 2882])

    # # get normalized bead coords
    # input_pts_frac = tfms.perform_coordinate_normalization(corners[nissl_id]['topleft_x'],
    #                                                        corners[nissl_id]['topleft_y'],
    #                                                        corners[nissl_id]['botright_x'],
    #                                                        corners[nissl_id]['botright_y'],
    #                                                        corners[nissl_id]['topright_x'],
    #                                                        corners[nissl_id]['topright_y'],
    #                                                        corners[nissl_id]['botleft_x'],
    #                                                        corners[nissl_id]['botleft_y'],
    #                                                        extents,
    #                                                        input_pts)

    # get normalized bead coords
    input_pts_frac = tfms.coordinate_normalization_in_padded_space(corners[nissl_id]['topleft_x'],
                                                                   corners[nissl_id]['topleft_y'],
                                                                   corners[nissl_id]['botright_x'],
                                                                   corners[nissl_id]['botright_y'],
                                                                   corners[nissl_id]['topright_x'],
                                                                   corners[nissl_id]['topright_y'],
                                                                   corners[nissl_id]['botleft_x'],
                                                                   corners[nissl_id]['botleft_y'],
                                                                   extents,
                                                                   input_pts)
    # input_pts_frac = tfms.perform_coordinate_normalization(0,
    #                                                        0,
    #                                                        0,
    #                                                        0,
    #                                                        0,
    #                                                        0,
    #                                                        0,
    #                                                        0,
    #                                                        extents,
    #                                                        input_pts)
    


    # get puck rescaled image coords

    # puck_rescaled_img_coords_map = map(lambda pt: [int(math.floor(pt[0]*img_width_ss_rescaled)), img_height_ss_rescaled-int(math.floor(pt[1]*img_height_ss_rescaled))], input_pts_frac)
    puck_rescaled_img_coords_map = map(lambda pt: [int(math.floor(pt[0]*img_width_ss_rescaled)), img_height_ss_rescaled-int(math.floor(pt[1]*img_height_ss_rescaled))], input_pts_frac)
    puck_rescaled_img_coords = list(puck_rescaled_img_coords_map)

    # get tfmed_ss
    PIX_SIZE = 0.000507
    puck_rescaled_world = list(map(lambda pt:[-pt[0]*PIX_SIZE, -pt[1]*PIX_SIZE], puck_rescaled_img_coords))

    # tfm_aff = tfms.get_histolozee_affine_tfm(hz_project_ss_file, nissl_id)
    tfm_aff = tfms.get_histolozee_affine_tfm_contructed(hz_project_ss_file, nissl_id)
    print(tfm_aff)

    # tfm_1 = tfms.get_affine_transform(0, [1,1], [0.438048, 0.438048])
    # tfm_1_inv = tfms.get_affine_transform(0, [1,1], [-0.438048, -0.438048])

    pts_tfmed_ss = []
    # transform to transformed ss world coordinates
    for pt in puck_rescaled_world:
            pt_homo = np.array([pt[0], pt[1], 1])
            pt_tfmed_ss = tfm_aff@pt_homo
            # print(pt_homo, pt_tfmed_ss)
            pts_tfmed_ss.append(pt_tfmed_ss[:-1])
            # pts_tfmed_ss.append(tfm_mat@pt)

    ## flipping x axis to account for nifty coordinates
    # pts_tfmed_ss = [ [pt[0], pt[1]] for pt in pts_tfmed_ss]


    # get tfmed_ss_normalized

    top = -0.0759877
    left = -0.0759877
    width = 1.02807
    height = 1.02807

    pts_tfmed_ss_normalized = [[(-pt[0]-left)/width, abs(pt[1]+top)/height] for pt in pts_tfmed_ss]
    # pts_tfmed_ss_normalized = [[-(pt[0]-left)/width, -(pt[1]-top)/height] for pt in pts_tfmed_ss]

    # get tfmed_ss_nifti

    height=img_width_ss_tfmed # dimensions of transformed ss images exported from histolozee
    width=img_height_ss_tfmed
    # pts_tfmed_ss_nifti = [[-pt[0]*height, -pt[1]*width] for pt in pts_tfmed_ss_normalized]
    pts_tfmed_ss_nifti = [[-pt[0]*width, -pt[1]*height] for pt in pts_tfmed_ss_normalized]


    # warped_coords
    nis_idx = str(nissl_id).zfill(3)
    intrim_pts_file = tfm_folder+"/"+nis_idx+"_intrim_pts.txt"
    with open(intrim_pts_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for row in pts_tfmed_ss_nifti:
            writer.writerow(row)

    from_fiducials_file = tfm_folder+"/"+str(int(nis_idx))+"_f.csv"
    to_fiducials_file = tfm_folder+"/"+str(int(nis_idx))+"_t.csv"
    # nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
    temp_output_csv_file = tfm_folder+"/"+nis_idx+"_tmp_output.csv"
    # call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
    subprocess.run(["./build/cmapper2", 
                    from_fiducials_file,
                    to_fiducials_file,
                    intrim_pts_file,
                    "None",
                    nis_idx,
                    temp_output_csv_file])


    ## read warped coords from temp file

    warped_coords = []
    with open(temp_output_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [row[0], row[1]]
            warped_coords.append(point)


    # frac_coords
    width = img_width_nis_tfmed
    height = img_height_nis_tfmed

    frac_coords = [ [-float(pt[0])/width, -float(pt[1])/height] for pt in warped_coords ]

    # TODO removal of hard coding
    width = 2.93578
    height = 2.58504
    top = 0.964691 # TODO figure out reason for sign inversion in histolozee export popup!!
    left = 0.670088

    # hz_nis_world_coords
    hz_nis_world_coords = [[left-(pt[0]*width), top-(pt[1]*height)] for pt in frac_coords]

    # print(hz_nis_world_coords)

    # convert to chuck_space_img_coords
    chuck_world_img_coords = [[-int(math.floor((pt[0]-left)*(img_width_nis_tfmed/width))), -int(math.floor((pt[1]-top)*(img_height_nis_tfmed/height)))] for pt in hz_nis_world_coords]

    # print(chuck_world_img_coords)

    # print("top", top, "left", left, "height", height, "width", width)
    op_file = op_folder+"/chuck_sp_img_coords_"+nis_idx+".csv"
    with open(op_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        # for idx, row in enumerate(pts_tfmed_ss):
        # for idx, row in enumerate(pts_tfmed_ss_normalized):
        # for idx, row in enumerate(puck_rescaled_img_coords):
        for idx, row in enumerate(chuck_world_img_coords):
            writer.writerow(row)


