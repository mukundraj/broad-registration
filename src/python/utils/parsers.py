"""
Module with functions to help parse Histolozee project xml.

Created by Mukund on 2022-03-04

"""

import xml.dom.minidom as parser
import numpy as np
import json
import lxml.etree as ET
import os

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
    mapper_to_id = {}
    seg_labels = root.find("segmentation-labels")
    for srno, label in enumerate(seg_labels.iter("label")):
        id = label.get("ID")
        color = label.get("color").replace("#", "")
        mapper[color] = srno+100
        mapper_to_id[srno+100] = id
        print (id, color, mapper[color])


    return mapper, mapper_to_id

""" Reads in histolozee project xml and nissl idx and returns stored
corresponding affine transform.

Created by Mukund on 2022-03-15
"""
def get_tfm_from_hz_xml(hz_project_file, nis_idx):

    tree = ET.parse(hz_project_file)
    root = tree.getroot()


    tfm_aff = None
    for slide in root.iter("slide"):
        # print (slide.tag, slide.attrib)
        img_relative_path = slide.get("image")
        basename = os.path.basename(img_relative_path)
        elm_idx = slide.get("image").split("_")[1]
        if (int(elm_idx)==int(nis_idx)):
            # print("elm_idx:", elm_idx)
            tfm = slide.find("transformation")
            tfm = str(tfm[5])
            tfm_aff = get_np_array_from_tfm_string(tfm)

    return tfm_aff

