"""
Script to convert+rescale vsi nissl VSI files to tiff/png

Read all files, convert and rescale, and save with same filename with new
extension. Depends on bfconvert commandline tool, in turn on java runtime (see
references).

Created by Mukund on 2022-02-16

Usage: python script.py path_to_input_folder path_to_output_folder opfileformat

Usage example: 

python src/python/scripts/convert_from_vsi.py \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/vsi \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/png_from_vsi \
png

References:
- https://docs.openmicroscopy.org/bio-formats/5.7.3/users/comlinetools/conversion.html
- https://java.com/en/download/

"""
from os import listdir
from os.path import isfile, join
import sys
import os
import subprocess

input_path = sys.argv[1]
output_path = sys.argv[2]
opfileformat = sys.argv[3]

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

