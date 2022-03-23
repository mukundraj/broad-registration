
import transforms3d as tfm
import math
import numpy as np
import src.python.utils.parsers as parsers
import lxml.etree as ET

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


"""
Read histolozee slide seq xml and return transformation in affine matrix format

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


"""
Performs coordinate normalization including removal of padding

Inputs: corners, input_pts

Outputs: normalized_coordinates

Created by Mukund on 2022-03-21
"""
def perform_coordinate_normalization(topleft_x, topleft_y,
                                     botright_x, botright_y,
                                     topright_x, topright_y,
                                     botleft_x, botleft_y,
                                     input_pts):

    s = np.array([[topleft_x, topleft_y, 1],
                  [topright_x, topright_y, 1],
                  [botleft_x, botleft_y, 1],
                  [botright_x, botright_y, 1]], dtype=float).T

    d = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float).T

    M_rec,resid,rank,sing = np.linalg.lstsq(s.T,d.T)
    M_rec = M_rec.T

    # format input
    N = len(input_pts)
    pts = np.array(input_pts)
    pts = np.hstack((pts, np.ones(N).reshape(N,1))).T

    # transform input
    tfmed_pts = M_rec@pts

    tfmed_pts = list(tfmed_pts.T)

    return tfmed_pts



