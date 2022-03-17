
import transforms3d as tfm
import math
import numpy as np

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

