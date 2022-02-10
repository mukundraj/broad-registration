import os
from os import listdir
from os.path import isfile, join
import shutil
import subprocess

parent = "/Users/mraj/Desktop/work/data/mouse_atlas/initial_converted_subset"

files_to_copy = []
path_to_copy_to = "/Users/mraj/Desktop/work/data/mouse_atlas/initial_converted_subset_/"

for i in range(1,91):
    path =parent+"/"+str(i)
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    files = [ fi for fi in onlyfiles if fi.endswith(".tif") ]
    if (len(files)==1):
        filename  = parent+'/'+str(i)+'/'+files[0]
        files_to_copy.append(filename)
        dst_filename = path_to_copy_to+files[0]
        base = os.path.splitext(dst_filename)[0]
        png_filename = base+".png"

        print (filename, png_filename)
        # shutil.copy(filename, dst_filename)
        subprocess.run(["convert", filename, png_filename])
