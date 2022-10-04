"""
Script to convert+rescale vsi nissl VSI files to tiff

Read all files, convert and rescale, and save with same filename with new
extension. Depends on bfconvert commandline tool, in turn on java runtime (see
references).

Created by Mukund on 2022-02-16

Usage: python script.py path_to_input_folder path_to_output_folder path_to_output_forslicer

Usage example: 

python src/python/scripts/convert_from_vsi.py \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_raw_data/vsi \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/tiff_from_vsi \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/tiff_from_vsi_forslicer

References:
- https://docs.openmicroscopy.org/bio-formats/5.7.3/users/comlinetools/conversion.html
- https://java.com/en/download/
- parallization: can be easily parallized - https://shuzhanfan.github.io/2017/12/parallel-processing-python-subprocess/

"""
from os import listdir
from os.path import isfile, join
import sys
import os
import subprocess

input_path = sys.argv[1]
output_path = sys.argv[2]
output_path_forslicer = sys.argv[3]
opfileformat = "tiff"

# read files
onlyfiles = [f for f in listdir(input_path) if isfile(join(input_path, f))]
files = [f for f in onlyfiles if os.path.splitext(f)[-1]==".vsi"]

## convert and write out
for file in files:
    fullpath = input_path+"/"+file
    print(fullpath)
    fullpath_op = output_path+"/nis_"+file
    pre, ext = os.path.splitext(fullpath_op)
    fullpath_op = pre+"."+opfileformat
    print(fullpath_op)
    # subprocess.run(["/Users/mraj/Downloads/bftools/bfconvert", "-tilex", "512", "-tiley", "512", fullpath , fullpath_op])
    subprocess.run(["/Users/mraj/Desktop/work/pkgs/bftools/bfconvert", "-series", "3", fullpath , fullpath_op])

# convert to tif format needed by histolozee
subprocess.run(["fd", "^.*tiff$" , "-x", "convert", "{}", "-define", "tiff:tile-geometry=256x256", "ptif:{.}.tif" ], cwd=output_path)

# delete intermediate tiff files produced after initial conversion from vsi
# subprocess.run(["fd", "tiff$", "-X", "rm"], cwd=output_path)

subprocess.run(["fd", "tiff$", "-x", "mv", "{}", output_path_forslicer], cwd=output_path)



