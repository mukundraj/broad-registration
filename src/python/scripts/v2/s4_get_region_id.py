"""
Script to read in list of puck positions and produce a list of Chuck region ids.

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


# file = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/hz-project-ss.zee"


nis_idx = sys.argv[1]
hz_project_ss_file = sys.argv[2]
tfm_folder = sys.argv[3]
hz_project_file = sys.argv[4]
labelmap_folder = sys.argv[5]
ss_rescaled_folder = sys.argv[6]
nissl_folder = sys.argv[7]
mapper_to_id_file = sys.argv[8]
input_csv_file = sys.argv[9]
output_csv_file = sys.argv[10]

# import xml.etree.ElementTree as ET
import lxml.etree as ET
tree = ET.parse(hz_project_ss_file)
root = tree.getroot()


tfm_aff = None
for slide in root.iter("slide"):
    # print (slide.tag, slide.attrib)
    elm_idx = slide.get("image").split("_")[4]
    if (elm_idx==nis_idx):
        print(elm_idx)
        tfm = slide.find("transformation")
        tfm = str(tfm[5])
        tfm_aff = parsers.get_np_array_from_tfm_string(tfm)

        # pt = np.array([1,1,1])

        # print(np.matmul(tfm_aff, pt))



# read input puck fractions
input_pts_frac = []
with open(input_csv_file, newline='\n') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # print(row['first_name'], row['last_name'])
        row = list(map(float, row))
        point = [row[0], row[1]]
        input_pts_frac.append(point)

# width = 1728
# height = 1728
# test_inputs = [[304,896], [316, 1157], [836, 1307], [319, 549]] # output 114 112 113 100
# input_pts_frac = [[pt[0]/width, pt[1]/height] for pt in test_inputs]

print("Input points as fractions:")
print(input_pts_frac)

# transoform to puck rescaled image coords
ss_file_name = ss_rescaled_folder+"/ss_"+nis_idx+"_hz.tif"
img = Image.open(ss_file_name)

print("Img width and height:")
print(img.width,img.height)

puck_rescaled_img_coords_map = map(lambda pt: [int(math.floor(pt[0]*img.width)), img.height-int(math.floor(pt[1]*img.height))], input_pts_frac)
puck_rescaled_img_coords = list(puck_rescaled_img_coords_map)


print("Puck rescaled ss img coords:")
print(puck_rescaled_img_coords, "\n")

# transform to puck rescaled world coords

PIX_SIZE = 0.000507
puck_rescaled_world = list(map(lambda pt:[-pt[0]*PIX_SIZE, -pt[1]*PIX_SIZE], puck_rescaled_img_coords))

print("Puck rescaled ss world coords:")
print(puck_rescaled_world, "\n")

pts_tfmed_ss = []


# print( tfm_aff, "\n")
# tfm_aff = tfms.get_affine_transform(angle_deg=90.12, scaling=[1,-1], translation=[0,0])
# print( tfm_aff, "\n")

tfm_mat = tfms.get_SR_matrix(angle_deg=90.12, scaling=[1,-1])
# print(tfm_mat)
# transform to transformed ss world coordinates
for pt in puck_rescaled_world:
        pt_homo = np.array([pt[0], pt[1], 1])
        pt_tfmed_ss = tfm_aff@pt_homo
        # print(pt_homo, pt_tfmed_ss)
        pts_tfmed_ss.append(pt_tfmed_ss[:-1])
        # pts_tfmed_ss.append(tfm_mat@pt)

## flipping x axis to account for nifty coordinates
pts_tfmed_ss = [ [pt[0], pt[1]] for pt in pts_tfmed_ss]
print("Transformed ss world coordinates:")
print(pts_tfmed_ss, "\n")

# transform to transformed ss normalized coordinates using **histozee bounds**
top = -0.0759877
left = -0.0759877
width = 1.02807
height = 1.02807

pts_tfmed_ss_normalized = [[(-pt[0]-left)/width, abs(pt[1]+top)/height] for pt in pts_tfmed_ss]

print("Transformed ss normalized coordinates:")
print(pts_tfmed_ss_normalized)
print("\n")

# transform to transformed slide seq NIfTI  - using image dims, pixel resolution, and multiply by minus
height=2048
width=2048
# pts_tfmed_ss_nifti = [[-pt[0]*height, -pt[1]*width] for pt in pts_tfmed_ss_normalized]
pts_tfmed_ss_nifti = [[-pt[0]*width, -pt[1]*height] for pt in pts_tfmed_ss_normalized]

print("Transformed slide seq nifti coords:")
print(pts_tfmed_ss_nifti)

# transform to warped world coords

intrim_pts_file = tfm_folder+"/"+nis_idx+"_intrim_pts.txt"
with open(intrim_pts_file, 'w', newline='\n') as csvfile:
    writer = csv.writer(csvfile)
    for row in pts_tfmed_ss_nifti:
        writer.writerow(row)

from_fiducials_file = tfm_folder+"/"+str(nis_idx)+"_f.csv"
to_fiducials_file = tfm_folder+"/"+str(nis_idx)+"_t.csv"
nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
temp_output_csv_file = tfm_folder+"/"+nis_idx+"_tmp_output.csv"
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
print("\nFrac coords in chuck space:")
print(frac_coords)

# TODO removal of hard coding
width = 2.93578
height = 2.58504
top = 0.964691 # TODO figure out reason for sign inversion in histolozee export popup!!
left = 0.670088

hz_nis_world_coords = [[left-(pt[0]*width), top-(pt[1]*height)] for pt in frac_coords]
# hz_nis_world_coords.append([-0.720163,  -0.321557]) # aligned
# hz_nis_world_coords.append([-0.64801, -1.14049]) # original
# hz_nis_world_coords.append([1,0]) # test
# hz_nis_world_coords.append([-1.34186, -0.206216]) # test
# hz_nis_world_coords.append([-1.04103, 0.235432]) # test
# hz_nis_world_coords.append([-1.18825, 0.267435]) # test

print("\nhz_nis_world_coords")
print(hz_nis_world_coords)

## get transform
tfm_aff = parsers.get_tfm_from_hz_xml(hz_project_file=hz_project_file, nis_idx=nis_idx)


# tfm_aff = tfms.get_affine_transform(90, [1,1], [0.824636, 1.177000])

# tfm_1 = tfms.get_affine_transform(0, [1,1], [1,0])
# tfm_2 = tfms.get_affine_transform(90, [1,1], [0, 0])
# tfm_3 = tfms.get_affine_transform( 0, [1,1], [-1, 0])

tfm_1 = tfms.get_affine_transform(0, [1,1], [0.824636, 1.177000])
tfm_2 = tfms.get_affine_transform(268.265015, [1,1], [0, 0])
tfm_3 = tfms.get_affine_transform( 0, [1,1],[-0.824636, -1.177000])
tfm_4 = tfms.get_affine_transform( 0, [1,1],[0.067142, 1.034830])

# tfm_1 = tfms.get_affine_transform( 0, [1,1],[-0.067142, -1.034830])
# tfm_2 = tfms.get_affine_transform(0, [1,1], [0.824636, 1.177000])
# tfm_3 = tfms.get_affine_transform(-268.265015, [1,1], [0, 0])
# tfm_4 = tfms.get_affine_transform( 0, [1,1],[-0.824636, -1.177000])

tfm_aff = tfm_4@tfm_3@tfm_2@tfm_1

## invert transform
tfm_aff = np.linalg.inv(tfm_aff)

# tfm_aff = tfm_2

## apply transform to go to original nissl world space
pts_nissl_world = []
for pt in hz_nis_world_coords:
        pt_homo = np.array([pt[0], pt[1], 1])
        pt_tfmed_ss = tfm_aff@pt_homo
        # print(pt_homo, pt_tfmed_ss)
        pts_nissl_world.append(pt_tfmed_ss[:-1])
        # pts_tfmed_ss.append(tfm_mat@pt)

print("\nOriginal nissl world coords:")
print(pts_nissl_world)

## convert to original nissl fractional coords from hz nissl world coords
# nissl_bb_left = -0.670088
# nissl_bb_top = -0.964691
# width = 2.93578
# height = 2.58504

nissl_bb_left = 0
nissl_bb_top = 0
width = 0.000507*3253
height = 0.000507*4643

# convert to fractional coordinates of original nissl
nissl_frac_coords = [[(-nissl_bb_left-pt[0])/width, (nissl_bb_top-pt[1])/height] for pt in pts_nissl_world]

print("\nNissl fractional coords:")
print(nissl_frac_coords)

# convert to image space of original nissl
##TODO fix hardcoding of height and width below
height = 4643
width = 3253 

nissl_img_corrds = [[int(math.floor(pt[0]*width)), int(math.floor(pt[1]*height))] for pt in nissl_frac_coords]

out_of_bounds = [idx for idx in range(len(nissl_img_corrds)) if nissl_img_corrds[idx][0]>=width or nissl_img_corrds[idx][1]>=height]



print("\nNissl img coords:")
print(nissl_img_corrds)

out_of_bounds_coords = []
for idx in out_of_bounds:
    out_of_bounds_coords.append([nissl_img_corrds[idx], nissl_img_corrds[idx]])
    nissl_img_corrds[idx] = [1,1]

print("out of bounds")
print(out_of_bounds_coords)
print(len(out_of_bounds))

# sample from nrrd
nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
readdata, header = nrrd.read(nrrd_path)
print(readdata.shape)
print(header)

labels = [readdata[pt[0]][pt[1]] for pt in nissl_img_corrds]

dbfile = open('output/mapper_to_id.pickle', 'rb')
mapper_to_id = pickle.load(dbfile)
mapper_to_id[0] = "NOLABEL"
mapper_to_id[-1] = "OUTDOM"
# write output file
with open(output_csv_file, 'w', newline='\n') as csvfile:
    writer = csv.writer(csvfile)
    for idx, row in enumerate(labels):
        if (idx in out_of_bounds):
            row = -1
        line = [row, mapper_to_id[row], nissl_img_corrds[idx][0], nissl_img_corrds[idx][1]]
        # line = [row, mapper_to_id[row]]
        writer.writerow(line)
        print(line)


# annotations=["00","10","01","11", "0.5,0"]
# plt.scatter([pt[0] for pt in pts_tfmed_ss_nifti], [pt[1] for pt in pts_tfmed_ss_nifti])
# for i, label in enumerate(annotations):
#     plt.annotate(label, (pts_tfmed_ss_normalized[i][0], pts_tfmed_ss_normalized[i][1]))
# plt.axis('equal')
# plt.show()

# doc = parser.parse(file)
# print (doc.nodeName)
# print (doc.firstChild.tagName)

# slides = doc.getElementsByTagName("slide")

# tfms = slides[0].getElementsByTagName("transformation")
# print(slides[0].getAttribute("image"))

# for tfm in tfms:
#     print ("rotation: ", tfm.getElementsByTagName("rotation")[0].getAttribute("k"))
#     comments = tfm.xpath('/comment()')
#     print (comments)
# print(len(slides))

# print("done")
