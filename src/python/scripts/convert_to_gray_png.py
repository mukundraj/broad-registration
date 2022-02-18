"""

Script to read in png/tiff files in a folder and convert to grayscale pngs

Usage example:

python ./src/python/scripts/convert_to_gray_png.py \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/tif_8x_downsampled \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/tif_8x_downsampled_gray

"""
import os
from os import listdir
from os.path import isfile, join
import shutil
import subprocess
import sys


# parent = "/Users/mraj/Desktop/work/data/mouse_atlas/initial_converted_subset"
# parent = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/tif_8x_downsampled"
parent = sys.argv[1]

files_to_copy = []
# path_to_copy_to = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/tif_8x_downsampled_gray"
path_to_copy_to = sys.argv[2]

# for i in range(1,91):
#     path =parent+"/"+str(i)
#     onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
#     files = [ fi for fi in onlyfiles if fi.endswith(".tif") ]
#     if (len(files)==1):
#         filename  = parent+'/'+str(i)+'/'+files[0]
#         files_to_copy.append(filename)
#         dst_filename = path_to_copy_to+files[0]
#         base = os.path.splitext(dst_filename)[0]
#         png_filename = base+".png"
#         print (filename, png_filename)
#         # shutil.copy(filename, dst_filename)
#         subprocess.run(["convert", filename, "-set", "colorspace", "Gray", "-separate", "-average", png_filename])

path = parent

onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
files = [ fi for fi in onlyfiles if fi.endswith(".tif") or fi.endswith(".png") ]
for i in range(len(files)):
    filename =parent+'/'+files[i]
    print (filename)
    base = os.path.splitext(filename)[0]
    # basefilename = base.replace("downsampled", "downsampled_gray")
    basefilename = os.path.basename(base)
    png_filename = path_to_copy_to+"/"+basefilename+".png"
    print(png_filename)
    # if (len(files)==1):
    #     filename  = parent+'/'+str(i)+'/'+files[0]
    #     files_to_copy.append(filename)
    #     dst_filename = path_to_copy_to+files[0]
    #     base = os.path.splitext(dst_filename)[0]
    #     png_filename = base+".png"
    #     print (filename, png_filename)
    #     # shutil.copy(filename, dst_filename)
    subprocess.run(["convert", filename, "-set", "colorspace", "Gray", "-separate", "-average", png_filename])
