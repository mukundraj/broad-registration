"""
Script to rescale slide seq files in folder to 30% smaller and sets the image
spacing to 1mm to match with Nissl slices. Reads paths from config file.

Usage : python src/python/scripts/v2/s0_resize_ss.py path_to_config_file

Usage example:

python src/python/scripts/v2/s0a_resize_ss.py config.yml

Created by Mukund on 2022-03-02
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

config_file = sys.argv[1]
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    conf = config["slide_seq"]

ip_path = conf["ip_path_ss"]
op_path = conf["op_rescaled_ss"]
# op_path_forslicer = conf["op_rescaled_ss_forslicer"]

path = ip_path

onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
files = [ fi for fi in onlyfiles if fi.endswith(".png") ]

print(files)

for i in range(len(files)):
    # adding prefix to adress issue of slicer reading entire sequence otherwise
    ip_filename = ip_path + '/' + files[i]
    op_filename = op_path + '/' + "ss_"+files[i]
    print (ip_filename)

    # op_filename = op_filename.replace("slide_seq_imgs", "slide_seq_imgs_rescaled")
    padded_id = os.path.basename(op_filename).split(".")[0].replace("ss_","").zfill(3)
    op_filename = op_path+"/ss_"+str(padded_id)+".tiff"
    print(op_filename)

    subprocess.run(["convert", ip_filename, "-resize", "30%", "-density", "72", "-units",  "pixelsperinch", op_filename])

    # Read into ITK file
    reader = sitk.ImageFileReader()
    # reader.SetImageIO("TIFFImageIO")
    reader.SetFileName(op_filename)
    image = reader.Execute()
    image.SetSpacing((1, 1))

    writer = sitk.ImageFileWriter()
    writer.SetFileName(op_filename)
    writer.Execute(image)

    # subprocess.run(["fd", "^.*tiff$" , "-x", "convert", "{}", "-define", "tiff:tile-geometry=256x256", "ptif:{.}.tif" ], cwd=output_path)



# copy files to forslicer folder
# subprocess.run(["fd", "-x", "cp", "{}", op_path_forslicer+"/"+"{}"], cwd=op_path)

# convert files in first folder to histolozee tif format
subprocess.run(["fd", "^.*tiff$" , "-x", "convert", "{}", "-define", "tiff:tile-geometry=256x256", "ptif:{.}_hz.tif" ], cwd=op_path)

# clear first folder
subprocess.run(["fd", "-e", "tiff", "-X", "rm"], cwd=op_path)

# copy hz tiffs to first folder
# subprocess.run(["fd", "^.*hz.tif$", "-x", "mv", "{}", op_path], cwd=op_path_forslicer)


