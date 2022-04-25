""" Script to read in list of puck positions and produce a list of Chuck region
ids. Use src/python/gists/render_tfmed_poses.py for redering images for
visually verifying output.

Usage:

python script.py\
    nis_idx (nissl slide index)\
    path to histolozee project xml with slide-seq initial transforms\
    path to folder with transform files ({nis_idx}_f.csv {nis_idx}_t.csv)\
    path to histolozee project xml with nissl initial transforms\
    path to labelmap nrrds folder (lmap{nis_idx}.nrrd)\
    path to slide seq scaled images\
    path to nissl original images\
    path to mapper to chuck segment ids \
    path to input csv with input bead coordinates as fraction of image width/height\
    path to output.csv - an output filename to store segment ids\

python s4_get_region_id.py \
    inp: nissl_id \
    inp: path to slide-seq histolozee .zee file
    inp: path to first set of slicer transform (nissl/ss alignment)
    inp: path to nissl histolozee .zee file
    inp: path to labelmap nrrds folder
    inp: path to folder with bead coords
    inp: path to file mapping old and new nissl names (needed to determine labelmap nrrd filenames)
    inp: path to corners file storing extent of bead position in ggpot figures
    inp inp: dimensions of rescaled slide seq image
    out: path to output folder to store csv files mapping beads to segment ids

Usage example:

python src/python/scripts/v2/s4_get_region_id.py \
    143 \
    "/Users/mraj/Desktop/sample-hz/hz-project-ss.zee" \
    "/Users/mraj/Desktop/sample-hz/transforms" \
    "/Users/mraj/Desktop/sample-hz/sample.zee" \
    "/Users/mraj/Desktop/sample-hz/labelmaps" \
    "/Users/mraj/Desktop/sample-hz/slide_seq_imgs_rescaled" \
    "/Users/mraj/Desktop/sample-hz"\
    "/Users/mraj/Desktop/sample-hz/mapper_to_id.pickle" \
    "/Users/mraj/Desktop/sample-hz/input.csv" \
    "/Users/mraj/Desktop/sample-hz/output.csv"

python src/python/scripts/v2/s4_get_region_id.py \
    -1 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/hz-project-ss.zee \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s3_registered_ss/transforms \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/hz-project.zee \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s2_seg_ids/seg_output \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_raw_data/bead_coords \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s2_seg_ids/filenames_map.csv \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/corners.csv \
    1728 1728 \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s4_bead_to_segid/bead_to_segid

Created by Mukund on 2022-03-11

"""
import matplotlib.pyplot as plt
from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
from PIL import Image

# import xml.dom.minidom as parser
import numpy as np
import math
# import json

import src.python.utils.parsers as parsers
import csv

import transforms3d as t3d
import src.python.utils.transforms as tfms
import subprocess
import nrrd
import pickle
import src.python.utils.io as io


# file = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/hz-project-ss.zee"


nissl_id = int(sys.argv[1])
hz_project_ss_file = sys.argv[2]
tfm_folder = sys.argv[3]
hz_project_file = sys.argv[4]
labelmap_folder = sys.argv[5]
bead_coords_folder = sys.argv[6]
mapperfile_csv = sys.argv[7]
corners_file = sys.argv[8]
img_width_ss_tfmed = int(sys.argv[9])
img_height_ss_tfmed = int(sys.argv[10])
output_csv_folder = sys.argv[11]

# ss_rescaled_folder = sys.argv[6]
# nissl_folder = sys.argv[7]
# mapper_to_id_file = sys.argv[8]

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

# import xml.etree.ElementTree as ET
import lxml.etree as ET
tree = ET.parse(hz_project_ss_file)
root = tree.getroot()

img_dims = io.get_img_dimensions(mapperfile_csv)
mapper, mapper_to_id = parsers.get_label_dict(hz_project_file)

# corners_file = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/corners.csv"
# read corners csv
corners = io.get_corners_info(corners_file)

for nissl_id in nissl_ids:

    # bead_coords_folder = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_raw_data/bead_coords"
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
    print(extents)
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

    nis_idx = str(nissl_id).zfill(3)


    # input_csv_file = f'{input_csv_folder}/normed_coords_{nis_idx}.csv'
    # read input puck fractions
    # input_pts_frac = []
    # with open(input_csv_file, newline='\n') as csvfile:
    #     reader = csv.reader(csvfile)
    #     for row in reader:
    #         # print(row['first_name'], row['last_name'])
    #         row = list(map(float, row))
    #         point = [row[0], row[1]]
    #         input_pts_frac.append(point)




    print(f'Input points as fractions loaded for {nis_idx}')

    puck_rescaled_img_coords_map = map(lambda pt: [int(math.floor(pt[0]*img_width_ss_tfmed)), img_height_ss_tfmed-int(math.floor(pt[1]*img_height_ss_tfmed))], input_pts_frac)
    puck_rescaled_img_coords = list(puck_rescaled_img_coords_map)


    print(f'Puck rescaled ss img coords computed for {nis_idx}')

    # transform to puck rescaled world coords

    PIX_SIZE = 0.000507
    puck_rescaled_world = list(map(lambda pt:[-pt[0]*PIX_SIZE, -pt[1]*PIX_SIZE], puck_rescaled_img_coords))

    print(f'Puck rescaled ss world coords computed for {nis_idx}')

    pts_tfmed_ss = []

    # print( tfm_aff, "\n")
    # print( tfm_aff, "\n")
    # tfm_aff = tfms.get_affine_transform(angle_deg=90.12, scaling=[1,-1], translation=[0,0])

    # tfm_mat = tfms.get_SR_matrix(angle_deg=90.12, scaling=[1,-1])

    # tfm_mat2 = tfms.get_histolozee_affine_tfm_contructed2(hz_project_ss_file, nissl_id)

    # tfm_aff = parsers.get_tfm_from_hz_ss_xml(hz_project_ss_file, nis_idx)
    # tval = str(67;).zfill(3)
    # print(tval)
    print("nis_idx ", nis_idx)
    tfm_aff = tfms.get_histolozee_affine_tfm_contructed2(hz_project_ss_file=hz_project_ss_file, nis_idx=nis_idx)
    print(tfm_aff)
    # transform to transformed ss world coordinates
    for pt in puck_rescaled_world:
            pt_homo = np.array([pt[0], pt[1], 1])
            pt_tfmed_ss = tfm_aff@pt_homo
            # print(pt_homo, pt_tfmed_ss)
            pts_tfmed_ss.append(pt_tfmed_ss[:-1])
            # pts_tfmed_ss.append(pt_homo[:-1])
            # pts_tfmed_ss.append(tfm_mat@pt)


    ## flipping x axis to account for nifty coordinates
    pts_tfmed_ss = [ [pt[0], pt[1]] for pt in pts_tfmed_ss]
    print(f'Transformed ss world coordinates for {nis_idx}')
    # print(pts_tfmed_ss, "\n")

    # transform to transformed ss normalized coordinates using **histozee bounds**
    top = -0.0759877
    left = -0.0759877
    width = 1.02807
    height = 1.02807

    # pts_tfmed_ss_normalized = [[(-pt[0]+left)/width, (-pt[1]+top)/height] for pt in pts_tfmed_ss]
    pts_tfmed_ss_normalized = [[(-pt[0]-left)/width, abs(pt[1]+top)/height] for pt in pts_tfmed_ss]

    print(np.amax(np.array(pts_tfmed_ss_normalized), axis=0), np.amin(np.array(pts_tfmed_ss_normalized), axis=0))

    print(f'Transformed ss normalized coordinates {nis_idx}')
    # print(pts_tfmed_ss_normalized)
    # print("\n")

    # transform to transformed slide seq NIfTI  - using image dims, pixel resolution, and multiply by minus
    height=2048
    width=2048
    # pts_tfmed_ss_nifti = [[-pt[0]*height, -pt[1]*width] for pt in pts_tfmed_ss_normalized]
    pts_tfmed_ss_nifti_pos = [[pt[0]*width, pt[1]*height] for pt in pts_tfmed_ss_normalized]
    pts_tfmed_ss_nifti = [[-pt[0]*width, -pt[1]*height] for pt in pts_tfmed_ss_normalized]

    print(f'Transformed slide seq nifti coords for {nis_idx}')
    # print(pts_tfmed_ss_nifti)

    # transform to warped world coords

    intrim_pts_file = tfm_folder+"/"+nis_idx+"_intrim_pts.txt"
    with open(intrim_pts_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for row in pts_tfmed_ss_nifti:
            writer.writerow(row)


    from_fiducials_file = f'{tfm_folder}/{str(nissl_id)}_f.csv'
    to_fiducials_file = f'{tfm_folder}/{str(nissl_id)}_t.csv'
    nrrd_path = "RANDOMSTRING"
    temp_output_csv_file = f'{tfm_folder}/{str(nissl_id)}_tmp_output.csv'
    # call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
    subprocess.run(["./build/cmapper2", 
                    from_fiducials_file,
                    to_fiducials_file,
                    intrim_pts_file,
                    nrrd_path,
                    nis_idx,
                    temp_output_csv_file])


    # transform to chuck world coordinates

    ## read warped coords from temp file

    warped_coords = []
    with open(temp_output_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [row[0], row[1]]
            warped_coords.append(point)

    # warped_coords.append([-2394, -1008])
    ## convert to nissl hz world coords
    width = 4096
    height = 3606

    frac_coords = [ [-float(pt[0])/width, -float(pt[1])/height] for pt in warped_coords ]
    warped_coords_pos = [ [-float(pt[0]), -float(pt[1])] for pt in warped_coords ] # OK
    print(f'Frac coords in chuck space computed for {nis_idx}')
    # print(frac_coords)

    # TODO removal of hard coding
    width = 2.93578
    height = 2.58504
    top = 0.964691 # TODO figure out reason for sign inversion in histolozee export popup!!
    left = 0.670088

    hz_nis_world_coords = [[left-(pt[0]*width), top-(pt[1]*height)] for pt in frac_coords]

    print(f'hz_nis_world_coords for {nis_idx}')
    # print(hz_nis_world_coords)

    print(hz_project_file)

    ## get transform
    # tfm_aff_init = parsers.get_tfm_from_hz_xml(hz_project_file=hz_project_file, nis_idx=str(nissl_id))

    # tfm_1 = tfms.get_affine_transform(0, [1,1], [0.824636, 1.177000])
    # tfm_2 = tfms.get_affine_transform(268.265015, [1,1], [0, 0])
    # tfm_3 = tfms.get_affine_transform( 0, [1,1],[-0.824636, -1.177000])
    # tfm_4 = tfms.get_affine_transform( 0, [1,1],[0.067142, 1.034830])

    # tfm_aff = tfm_4@tfm_3@tfm_2@tfm_1

    # print(tfm_aff_init)

    # print(tfm_aff)
    mapper_to_new_filename, mapper_to_old_filename = io.get_filenames_map(mapperfile_csv)
    old_nis_filename = mapper_to_old_filename[nissl_id]
    old_nis_id = old_nis_filename.split("_")
    if(len(old_nis_id)==3):
        old_nis_id = int(old_nis_id[1])
    elif (old_nis_id==2):
        old_nis_id = int(old_nis_id.split(".")[0])
    else:
        assert(False)
    print("old_nis_id: ", old_nis_id)
    tfm_aff = tfms.get_histolozee_affine_tfm_contructed2(hz_project_ss_file=hz_project_file, nis_idx=str(old_nis_id))
    ## invert transform
    tfm_aff = np.linalg.inv(np.array(tfm_aff))

    # tfm_aff = tfm_2

    ## apply transform to go to original nissl world space
    pts_nissl_world = []
    for pt in hz_nis_world_coords:
            pt_homo = np.array([pt[0], pt[1], 1])
            pt_tfmed_ss = tfm_aff@pt_homo
            # print(pt_homo, pt_tfmed_ss)
            pts_nissl_world.append(pt_tfmed_ss[:-1])
            # pts_tfmed_ss.append(tfm_mat@pt)

    print(f'\nOriginal nissl world coords computed for ')
    # print(pts_nissl_world)


    ## convert to original nissl fractional coords from hz nissl world coords
    # nissl_bb_left = -0.670088
    # nissl_bb_top = -0.964691
    # width = 2.93578
    # height = 2.58504

    nissl_bb_left = 0
    nissl_bb_top = 0
    # width = 0.000507*3253
    # height = 0.000507*4643
    width = 0.000507*img_dims[nissl_id][0]
    height = 0.000507*img_dims[nissl_id][1]

    print("Height Width", height, width)
    # convert to fractional coordinates of original nissl
    nissl_frac_coords = [[(-nissl_bb_left-pt[0])/width, (nissl_bb_top-pt[1])/height] for pt in pts_nissl_world]

    print(f'\nNissl fractional coords computed for {nis_idx}')


    # convert to image space of original nissl
    ##TODO fix hardcoding of height and width below
    # height = 4643
    # width = 3253

    height = img_dims[nissl_id][1]
    width = img_dims[nissl_id][0]
    print("height width", height, width)

    nissl_img_corrds = [[int(math.floor(pt[0]*width)), int(math.floor(pt[1]*height))] for pt in nissl_frac_coords]

    out_of_bounds = [idx for idx in range(len(nissl_img_corrds)) if nissl_img_corrds[idx][0]>=width or nissl_img_corrds[idx][1]>=height]

    print(f'\nNissl img coords {nis_idx}')
    # print(nissl_img_corrds)

    out_of_bounds_coords = []
    for idx in out_of_bounds:
        out_of_bounds_coords.append([nissl_img_corrds[idx], nissl_img_corrds[idx]])
        nissl_img_corrds[idx] = [1,1]

    print(f'out of bounds')
    # print(out_of_bounds_coords)
    print(f'Len out of bounds {len(out_of_bounds)}')

    # sample from nrrd
    nrrd_filename = old_nis_filename.replace("tif", "nrrd")
    print(nrrd_filename)

    nrrd_path = f'{labelmap_folder}/{nrrd_filename}'
    readdata, header = nrrd.read(nrrd_path)
    print(readdata.shape)
    print(header)

    labels = [readdata[pt[0]][pt[1]] for pt in nissl_img_corrds]

    # dbfile = open('output/mapper_to_id.pickle', 'rb')
    # mapper_to_id = pickle.load(dbfile)
    mapper_to_id[0] = "NOLABEL"
    mapper_to_id[-1] = "OUTDOM"

    output_csv_file = f'{output_csv_folder}/bead_to_segid_{nis_idx}.csv'
    # write output file
    with open(output_csv_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx, row in enumerate(labels):
            if (idx in out_of_bounds):
                row = -1
            # line = [row, mapper_to_id[row], puck_rescaled_img_coords[idx][0], puck_rescaled_img_coords[idx][1]] # OK
            # line = [row, mapper_to_id[row], pts_tfmed_ss_nifti_pos[idx][0], pts_tfmed_ss_nifti_pos[idx][1]] # OK
            # line = [row, mapper_to_id[row], warped_coords_pos[idx][0], warped_coords_pos[idx][1]] # OK
            line = [row, mapper_to_id[row], nissl_img_corrds[idx][0], nissl_img_corrds[idx][1]]
            writer.writerow(line)


