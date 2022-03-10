""" Create a label map NRRD file from histolozee segmentation. Input are
labelmap tif images exported by Histolozee and Histolozee project xml file.
Label color values in project xml file need to be corrected since they don't
match with the output tifs.

Created by Mukund on 2022-02-28

Usage example:

python src/python/scripts/v2/s1_hz_seg_to_nrrd.py ./input/labels_config.yml group2

References:
- https://www.geeksforgeeks.org/how-to-convert-images-to-numpy-array/
- https://stackoverflow.com/questions/16414559/how-to-use-hex-without-0x-in-python
- https://stackoverflow.com/questions/42594695/how-to-apply-a-function-map-values-of-each-element-in-a-2d-numpy-array-matrix

"""
import sys
sys.path.append("src/python/utils")

import yaml
from PIL import Image
from numpy import asarray
import numpy as np
import os
import subprocess
import nrrd
# import "src/python/utils/parsers.py" as parsers
import parsers as parsers
# from src.python.utils import parsers

config_yaml = sys.argv[1]
group_id = sys.argv[2]

# read all filenames in group
with open(config_yaml) as file:
    files = yaml.load(file,  Loader=yaml.FullLoader)[group_id]
with open(config_yaml) as file:
    hz_xml_path = yaml.load(file, Loader=yaml.FullLoader)['hz-xml-'+group_id]
    print("\nxml-path:", hz_xml_path)

print("")
mapper = parsers.get_label_dict(hz_xml_path)
mapper['ffffff'] = 0

def mapper_wrapper(e,f,g ):
    # key = format(e,'x')+str(f)+str(g)
    key = format(e,'x')+format(f, 'x')+format(g, 'x')
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
    tiff_name = dirname+"/"+base.replace(".svg", "_labels.tif" )

    # print(tiff_name)
    # convert to tiff
    subprocess.run(["convert", file, tiff_name])

    # open tiff file
    img = Image.open(tiff_name)

    # convert img to numpy array
    nd = asarray(img)
    nd = np.swapaxes(nd,0,1)

    labels_current = lfunc(nd[:,:,0], nd[:,:,1], nd[:,:,2])
    ctr = ctr+1

    if (numpydata is not None):
        numpydata = np.stack([numpydata,labels_current], axis=2)
    else:
        numpydata = labels_current

print("")
# create nrrd file
nrrdfile = "output/labelmap.nrrd"
print ("Output labelmap written to: ", nrrdfile)
nrrd.write(nrrdfile, numpydata)
