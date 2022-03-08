"""
Module with functions to help parse Histolozee project xml.

Created by Mukund on 2022-03-04

"""

import xml.dom.minidom as parser
import numpy as np
import json
import lxml.etree as ET

"""

Created by Mukund on 2022-03-04
"""
def get_np_array_from_tfm_string(tfm_aff):

    tfm_aff = tfm_aff.replace("<!--2D affine matrix: ", "")
    tfm_aff = tfm_aff.replace("-->", "")
    tfm_aff = tfm_aff.split(";")
    tfm_aff[0] = tfm_aff[0].replace("[", "")
    tfm_aff[2] = tfm_aff[2].replace("]", "")
    tfm_aff[0] = tfm_aff[0].replace(" ", ",")
    tfm_aff[1] = tfm_aff[1][1:].replace(" ", ",")
    tfm_aff[2] = tfm_aff[2][1:].replace(" ", ",")
    # tfm_aff = tfm_aff.replace(" ", ",")
    tfm_json = "[["+tfm_aff[0]+"],["+tfm_aff[1]+"],["+tfm_aff[2]+"]]"
    # print (tfm_json)
    tfm_aff = np.array(json.loads(tfm_json))
    return tfm_aff

""" Reads in histolozee project xml path and returns a dictionary with label
colors mapped to integers

Created by Mukund on 2022-03-07
"""

def get_label_dict(file):
    tree = ET.parse(file)
    root = tree.getroot()

    mapper = {}
    seg_labels = root.find("segmentation-labels")
    for srno, label in enumerate(seg_labels.iter("label")):
        id = label.get("ID")
        color = label.get("color").replace("#", "")
        mapper[color] = srno+100
        print (id, color, mapper[color])


    return mapper
