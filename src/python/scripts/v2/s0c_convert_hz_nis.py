"""
Transforms aligned histolozee nissls to slicer format. Also removes alpha channel.

Usage Example:

python src/python/scripts/v2/s0c_convert_hz_nis.py config.yml

Created by Mukund on 2022-03-09

"""

import os
import yaml
import sys
import subprocess
from os import listdir
from os.path import isfile, join
import shutil
import subprocess
import SimpleITK as sitk
import numpy as np
import csv

config_file = sys.argv[1]
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    conf = config["nissl"]

ip_path = conf["transformed_hz_tif"]
op_path = conf["transformed_slicer_tif"]


onlyfiles = [f for f in listdir(ip_path) if isfile(join(ip_path, f))]
files = [ fi for fi in onlyfiles if fi.endswith(".tif") ]


for i in range(len(files)):
    # adding prefix to adress issue of slicer reading entire sequence otherwise
    sr_id = files[i].replace("ss_","").replace(".tif", "")
    ip_filename = ip_path + '/' + files[i]
    op_filename = op_path + '/' + sr_id+"_sl.tif"
    print (ip_filename, op_filename)
    reader = sitk.ImageFileReader()
    # reader.SetImageIO("TIFFImageIO")
    reader.SetFileName(ip_filename)
    image = reader.Execute()
    image.SetSpacing((1, 1))
    print (op_filename)
    writer = sitk.ImageFileWriter()
    writer.SetFileName(op_filename)
    writer.Execute(image)

subprocess.run(["fd", "^.*tif$" , "-x", "convert", "{}", "-alpha", "off", "{}" ], cwd=op_path)
