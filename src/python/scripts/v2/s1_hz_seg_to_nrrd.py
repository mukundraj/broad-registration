""" Create a label map NRRD file from histolozee segmentation. Input are
labelmap tif images exported by Histolozee and Histolozee project xml file.
Label color values in project xml file need to be corrected since they don't
match with the output tifs.

Created by Mukund on 2022-02-28

Usage:

python filename.py yml_config_file histolozee_config_tag label_files_tag

Usage example:

python src/python/scripts/v2/s1_hz_seg_to_nrrd.py ./config.yml histolozee_xml labelsvgs 

References:
- https://www.geeksforgeeks.org/how-to-convert-images-to-numpy-array/
- https://stackoverflow.com/questions/16414559/how-to-use-hex-without-0x-in-python
- https://stackoverflow.com/questions/42594695/how-to-apply-a-function-map-values-of-each-element-in-a-2d-numpy-array-matrix

"""
from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import yaml
from PIL import Image
from numpy import asarray
import numpy as np
import os
import subprocess
import nrrd
# import "src/python/utils/parsers.py" as parsers
# import parsers as parsers
import src.python.utils.parsers as parsers
# from src.python.utils import parsers
import pickle

config_yaml = sys.argv[1]
histolozee_config_tag = sys.argv[2]
label_files_tag = sys.argv[3]

print(histolozee_config_tag)

# read all filenames in group
with open(config_yaml) as file:
    files = yaml.load(file,  Loader=yaml.FullLoader)[label_files_tag]
with open(config_yaml) as file:
    hz_xml_path = yaml.load(file, Loader=yaml.FullLoader)[histolozee_config_tag]
    print("\nxml-path:", hz_xml_path)

print("")
mapper, mapper_to_id = parsers.get_label_dict(hz_xml_path)
mapper['ffffff'] = 0

dbfile = open('output/mapper_to_id.pickle', 'ab')
pickle.dump(mapper_to_id, dbfile)
dbfile.close()

def mapper_wrapper(e,f,g ):
    # key = format(e,'x')+str(f)+str(g)

    # converting to hexadecimal values
    d1 = format(e,'x')
    d2 = format(f,'x')
    d3 = format(g,'x')

    # dealing with single digit entries
    d1 = d1.zfill(2)
    d2 = d2.zfill(2)
    d3 = d3.zfill(2)
    key = d1+d2+d3
    if (key in mapper):
        label = mapper[key]
    else:
        label = 0
    return (label)

print("")
lfunc = np.vectorize(mapper_wrapper)
numpydata = None

ctr = 0
for file in files:
    print ("Reading label file: ", file)
    # create tiff filename
    base = os.path.basename(file)
    dirname = os.path.dirname(file)
    png_name = dirname+"/"+base.replace(".svg", ".png" )
    # tif_name = dirname+"/"+base.replace(".svg", "_test.tiff" )

    print(png_name)
    # convert to tiff/png
    # subprocess.run(["convert", file, tiff_name])
    # uses inkscape since this conversion not handled correctly by imagemagick
    subprocess.run(["inkscape", "--export-type=png", file ])

    # open tiff file
    img = Image.open(png_name)
    # convert img to numpy array
    # nd = asarray(img)
    nd = np.asarray(img)
    nd = np.swapaxes(nd,0,1)

    labels_current = lfunc(nd[:,:,0], nd[:,:,1], nd[:,:,2])
    ctr = ctr+1

    # print("labels_current", np.shape(labels_current))
    print(labels_current)
    labels_current = np.asarray(labels_current)

    print("labels_current", np.shape(labels_current))

    # img = Image.fromarray(labels_current)
    img = Image.fromarray(np.uint8(labels_current) , 'L')
    # img.save(tif_name)

    numpydata = labels_current
    print(np.shape(np.where(numpydata==109)))
    print(np.shape(numpydata))

    # if (numpydata is not None):
    #     numpydata = np.stack([numpydata,labels_current], axis=2)
    # else:
    #     numpydata = labels_current

    print("")
    # create nrrd file
    nrrdfile = dirname+"/labelmaps/"+base.replace(".svg", ".nrrd").replace("nis", "lmap")
    print ("Output labelmap written to: ", nrrdfile)
    nrrd.write(nrrdfile, numpydata)


