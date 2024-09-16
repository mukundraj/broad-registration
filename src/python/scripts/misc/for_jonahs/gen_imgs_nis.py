"""

python src/python/scripts/misc/for_jonahs/gen_imgs_nis.py


"""

import sys
from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import os
import SimpleITK as sitk
import subprocess

import src.python.utils.filename_maps as fm

ipdir = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png_ofix"
opdir = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/cs_nissls"
opdir_tifs = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/tif_nissls"
opdir_spfix = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/nissls_spfix"

pidToSrno, _ = fm.get_pidToSrno()

pids = list(range(1,208,2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

for pid in pids:


    pid_padded = str(pid).zfill(3)
    print (pid, pidToSrno[ pid ], pid_padded)

    from_file = f'{ipdir}/nis_{pid_padded}'
    to_file = f'{opdir}/nis_{pidToSrno[pid]}'
    to_file_tif = f'{opdir_tifs}/nis_{pidToSrno[pid]}'
    to_file_spfix = f'{opdir_spfix}/nis_{pidToSrno[pid]}'

    # copy from to to
    print (f'cp {from_file}.png {to_file}.png')

    os.system(f'cp {from_file}.png {to_file}.png')

    # convert to tif
    os.system(f'magick {to_file}.png {to_file_tif}.tif')

    # make spacing 1 
    reader = sitk.ImageFileReader()
    # reader.SetImageIO("TIFFImageIO")
    reader.SetFileName(f'{to_file_tif}.tif')
    image = reader.Execute()
    image.SetSpacing((1, 1))
    writer = sitk.ImageFileWriter()
    writer.SetFileName(f'{to_file_spfix}.tif')
    writer.Execute(image)

    # turn off alpha

    subprocess.run(["fd", "^.*tif$" , "-x", "convert", "{}", "-alpha", "off", "{}" ], cwd=opdir_spfix)

