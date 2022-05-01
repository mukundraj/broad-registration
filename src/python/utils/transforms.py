
import transforms3d as tfm
import math
import numpy as np
import src.python.utils.parsers as parsers
import lxml.etree as ET
import os
import csv
import subprocess
import json

def chuck_sp_to_allen3d(in_pts, nissl_id, tfm_folder, qnii_json_file, tmp_storage):
    """
    Reads in points in chuck space and returns points transformed to Allen CCF
    space.

    """
    nis_id_str = str(nissl_id).zfill(3)
    # write back negative
    intrim_pts_file = f'{tmp_storage}/chuck_space_coords_neg_{nis_id_str}.csv'
    with open(intrim_pts_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for row in in_pts:
            # print (row)
            writer.writerow([-row[0], -row[1]])

    from_fiducials_file = tfm_folder+"/"+str(nissl_id)+"_f.csv"
    to_fiducials_file = tfm_folder+"/"+str(nissl_id)+"_t.csv"
    # nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
    print(tmp_storage)
    temp_output_csv_file = tmp_storage+"/qnii_coords_"+nis_id_str+".csv"
    # call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
    subprocess.run(["./build/cmapper2",
                    from_fiducials_file,
                    to_fiducials_file,
                    intrim_pts_file,
                    "None",
                    nis_id_str,
                    temp_output_csv_file])

    warped_coords = []
    with open(temp_output_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(float, row))
            point = [-row[0], -row[1]]
            warped_coords.append(point)

    # moving on to qnii transformation
    f = open(qnii_json_file)

    # returns JSON object as
    # a dictionary
    qnii_data_unsorted = json.load(f)["slices"]
    qnii_data = {}
    img_dims = {}

    for item in qnii_data_unsorted:
        nissl_id_inloop = int(item["filename"][4:7])
        anch = item["anchoring"]
        # print(nissl_id)
        if nissl_id_inloop > 208:
            break
        qnii_data[nissl_id_inloop] = np.array([[anch[3], anch[4], anch[5]],
                                        [anch[6], anch[7], anch[8]],
                                        [anch[0], anch[1], anch[2]]])
        img_dims[nissl_id_inloop] = {"height": item["height"], "width":item["width"]}

    img_width = 4096 # fixme - hardcoding
    img_width = 3096

    input_pts = warped_coords
    # convert from left/right coordinate from left of image to left of animal
    # input_pts = [[img_width-pt[0], pt[1]] for pt in input_pts] # checkme - why this is not needed here but needed in other parts e.g mapper/get get id code

    pts = np.array(input_pts)

    # create normalized, homogeneous coords
    # pts = pts/np.array
    print ("width: ", img_dims[nissl_id]["width"], "height:", img_dims[nissl_id]["height"])
    den = np.array([img_dims[nissl_id]["width"], img_dims[nissl_id]["height"], 1]).reshape((3,1))
    N = len(pts)
    ones = np.ones(N).reshape((N,1))
    pts = np.concatenate((pts, ones), axis=1)
    # print(np.shape(ones),np.shape(pts))
    # print("den", den)
    pts = pts.T/den
    pts = pts.T
    # print("Pts normalized, homogeous:\n")
    # print(pts)
    # print("\n")
    # print(np.amax(pts, axis=0), np.amin(pts, axis=0))

    pts = pts@qnii_data[nissl_id]
    # print("\nPts in physical space:\n")
    # print(pts)
    pts_physical = pts
    print("\n")
    print(np.amax(pts, axis=0), np.amin(pts, axis=0))

    # transform to allen space
    ones = np.ones(N).reshape((N,1))
    pts = np.concatenate((pts, ones), axis=1)

    # print("\nPts in physical space homogenized:\n")
    # print(pts)
    T_allen = np.array([[0, 0, 25, 0],
                        [-25, 0, 0, 0],
                        [0, -25, 0, 0],
                        [13175, 7975, 0, 1]])

    # to convert to Allen image space
    T_allen2 = np.array([[0.04, 0, 0, 0],
                        [0, 0.04, 0, 0],
                        [0, 0, 0.04, 0],
                        [0, 0, 0, 1]])

    pts_allen = pts@T_allen@T_allen2

    print("\n Allen coords max and min")
    print(np.amax(pts_allen, axis=0), np.amin(pts_allen, axis=0))
    print("")

    pts_allen = list(pts_allen)
    return pts_allen


def get_affine_transform(angle_deg, scaling, translation):


    T = [translation[0], translation[1]]
    theta = math.radians(angle_deg)
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    Z = [scaling[0], scaling[1]]

    return (tfm.affines.compose(T, R, Z))

def get_SR_matrix(angle_deg, scaling):

    theta = math.radians(angle_deg)
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])

    S = np.array([[scaling[0], 0],[0, scaling[1]]])

    RS = R@S
    return RS

"""
Rotation followed by translation affine matrix that operates on homologous coordinates.
"""
def get_TR_affine(angle_deg, translation):

    T = [translation[0], translation[1]]
    theta = math.radians(angle_deg)
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    Z = [-1, 1]

    TR = tfm.affines.compose(T, R, Z)
    return TR


""" Read histolozee slide seq xml and return transformation in affine matrix
format stored in xml.

Created by Mukund on 2022-03-22
"""

def get_histolozee_affine_tfm(hz_project_ss_file, nis_idx):

    tree = ET.parse(hz_project_ss_file)
    root = tree.getroot()


    tfm_aff = None

    for slide in root.iter("slide"):
        # print (slide.tag, slide.attrib)
        elm_idx = int(slide.get("image").split("_")[4])
        if (elm_idx==nis_idx):
            print("Reading hz ss tfm for niss_id: ", elm_idx)
            tfm = slide.find("transformation")
            tfm = str(tfm[5])
            tfm_aff = parsers.get_np_array_from_tfm_string(tfm)

    return tfm_aff


""" Read histolozee slide seq xml and return transformation matrix in affine
format by manually constructing it from scale and reflection componencts in
xml.

Created by Mukund on 2022-03-23

"""

def get_histolozee_affine_tfm_contructed(hz_project_ss_file, nis_idx):

    tree = ET.parse(hz_project_ss_file)
    root = tree.getroot()
    tfm_aff = None
    for slide in root.iter("slide"):
        # print (slide.tag, slide.attrib)
        elm_idx = int(slide.get("image").split("_")[4])
        if (elm_idx==nis_idx):
            tfm = slide.find("transformation")
            rot_deg = float(tfm[0].get("k"))
            scale_i = float(tfm[1].get("i"))
            scale_j = float(tfm[1].get("j"))
            center_i = float(tfm[4].get("i"))
            center_j = float(tfm[4].get("j"))

            print(f"Transforms for id {elm_idx}: rot {rot_deg}, scale_i {scale_i}, scale_j {scale_j}, center_i {center_i}, center_j {center_j}")
            t1 = get_affine_transform(0, [1, 1], [center_i,center_j])
            t2 = get_affine_transform(0, [scale_i, scale_j], [0, 0])
            t3 = get_affine_transform(rot_deg, [1, 1], [0, 0])
            t4 = get_affine_transform(0, [1, 1], [-center_i, -center_j])
            tfm_aff = t4@t3@t2@t1

    return tfm_aff

def get_histolozee_affine_tfm_contructed2(hz_project_ss_file, nis_idx):

    tree = ET.parse(hz_project_ss_file)
    root = tree.getroot()
    tfm_aff = None
    for slide in root.iter("slide"):
        # print (slide.tag, slide.attrib)
        # print(slide.get("image"))
        # elm_idx = int(slide.get("image").split("_")[4])
        img_relative_path = slide.get("image")
        basename = os.path.basename(img_relative_path)
        # elm_idx = slide.get("image").split("_")[1]
        elm_idx = basename.split("_")[1] # invariant to dirname path
        # print(elm_idx)
        if (elm_idx.find(".tif")>0):    # to deal with nissl names witout 01 suffix
            elm_idx = elm_idx.split(".")[0]
        if (elm_idx==nis_idx):
            print(elm_idx, nis_idx)
            tfm = slide.find("transformation")
            rot_deg = float(tfm[0].get("k"))
            scale_i = float(tfm[1].get("i"))
            scale_j = float(tfm[1].get("j"))
            translation_i = float(tfm[3].get("i"))
            translation_j = float(tfm[3].get("j"))
            center_i = float(tfm[4].get("i"))
            center_j = float(tfm[4].get("j"))

            print(f"Transforms for id {elm_idx}: rot {rot_deg}, scale_i {scale_i}, scale_j {scale_j}, center_i {center_i}, center_j {center_j} translation {translation_i} {translation_j}")
            t1 = get_affine_transform(0, [1, 1], [center_i,center_j])
            t2 = get_affine_transform(0, [scale_i, scale_j], [0, 0])
            t3 = get_affine_transform(rot_deg, [1, 1], [0, 0])
            t4 = get_affine_transform(0, [1, 1], [-center_i, -center_j])
            t5 = get_affine_transform(0, [1, 1], [-translation_i, -translation_j])

            tfm_aff = t5@t4@t3@t2@t1

    return tfm_aff

"""
Performs coordinate normalization including removal of padding

Inputs: corners, input_pts

Outputs: normalized_coordinates

Created by Mukund on 2022-03-21
"""

# def perform_coordinate_normalization(topleft_x, topleft_y,
#                                      botright_x, botright_y,
#                                      topright_x, topright_y,
#                                      botleft_x, botleft_y,
#                                      extents,
#                                      input_pts):

def perform_coordinate_normalization(topleft_x, topleft_y,
                                     botright_x, botright_y,
                                     topright_x, topright_y,
                                     botleft_x, botleft_y,
                                     extents,
                                     input_pts):

    # s = np.array([[topleft_x, topleft_y, 1],
    #               [topright_x, topright_y, 1],
    #               [botleft_x, botleft_y, 1],
    #               [botright_x, botright_y, 1]], dtype=float).T
    s = np.array([[extents['min_x'], extents['min_x'], 1],
                  [extents['max_x'], extents['min_y'], 1],
                  [extents['min_x'], extents['max_y'], 1],
                  [extents['max_x'], extents['max_y'], 1]], dtype=float).T

    # d = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float).T
    d = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float).T

    M_rec,resid,rank,sing = np.linalg.lstsq(s.T,d.T)
    M_rec = M_rec.T

    # format input
    N = len(input_pts)
    pts = np.array(input_pts)
    pts = np.hstack((pts, np.ones(N).reshape(N,1))).T

    # transform input
    tfmed_pts = M_rec@pts

    tfmed_pts = tfmed_pts[:-1]

    tfmed_pts = list(tfmed_pts.T)
    # tpt = np.array([268, 274, 1])
    # tpt = np.array([5409, 5497, 1])
    # print(M_rec@tpt)
    # exit(0)

    return tfmed_pts


def coordinate_normalization_in_padded_space(topleft_x, topleft_y,
                                             botright_x, botright_y,
                                             topright_x, topright_y,
                                             botleft_x, botleft_y,
                                             extents,
                                             input_pts):

    s = np.array([[extents['min_x'], extents['min_y'], 1],
                  [extents['max_x'], extents['min_y'], 1],
                  [extents['min_x'], extents['max_y'], 1],
                  [extents['max_x'], extents['max_y'], 1]], dtype=float).T

    # d = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float).T
    # d = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float).T
    # topleft, top right, bot left, bot right
    side = 5760
    # side = 1728 
    d = np.array([[float(topleft_x)/side, float(topleft_y)/side, 1], 
                  [float(topright_x)/side, float(topright_y)/side, 1],
                  [float(botleft_x)/side, float(botleft_y)/side, 1],
                  [float(botright_x)/side, float(botright_y)/side, 1]], dtype=float).T

    M_rec,resid,rank,sing = np.linalg.lstsq(s.T,d.T)
    M_rec = M_rec.T

    # format input
    N = len(input_pts)
    pts = np.array(input_pts)
    pts = np.hstack((pts, np.ones(N).reshape(N,1))).T

    # transform input
    tfmed_pts = M_rec@pts

    tfmed_pts = tfmed_pts[:-1]

    tfmed_pts = list(tfmed_pts.T)
    # tpt = np.array([268, 274, 1])
    # tpt = np.array([5409, 5497, 1])
    # print(M_rec@tpt)
    # exit(0)

    return tfmed_pts


