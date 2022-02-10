""" 
Script to rescale slide seq files in folder to 50% smaller. Reads paths
from config file.

Usage : python scr_resize_ss.py path_to_config_file

Created by Mukund on 2022-02-10
"""
import os
import yaml
import sys
import subprocess
from os import listdir
from os.path import isfile, join
import shutil
import subprocess

print (sys.argv[1])

config_file = sys.argv[1]
with open(config_file, "r") as f:
    config = yaml.safe_load(f)
    conf = config["slide_seq"]

ip_path = conf["ip_path_ss"]
op_path = conf["op_rescaled_ss"]




path = ip_path

onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
print (onlyfiles)
files = [ fi for fi in onlyfiles if fi.endswith(".png") ]

print(files)

for i in range(len(files)):
    # adding prefix to adress issue of slicer reading entire sequence otherwise
    ip_filename = ip_path + '/' + files[i]
    op_filename = ip_path + '/' + "ss_"+files[i]
    print (ip_filename)

    op_filename = op_filename.replace("slide_seq_imgs", "slide_seq_imgs_rescaled")
    print(op_filename)


    subprocess.run(["convert", ip_filename, "-resize", "50%", op_filename])
