"""
Use VisuAlign fiducials to convert transform image and coordinates

Usage:

python s6a_nonlinear.py \
    nissl ids to process \
    path to VisuAlign json \
    path to input images \
    path to bead coords to be transformed \
    path to slicer fiducials folder \
    path for intrim fiducial files converted to csv format \
    path to output images \
    path to output transformed bead coords

Usage example:

python src/python/scripts/v2/s6a_nonlinear.py \
    141 \
    /Users/mraj/Desktop/transformed_hz_png/March23_Full_CCF.json \
    /Users/mraj/Desktop/transformed_hz_png \
    /Users/mraj/Desktop/forgif/chuck_space_img_coords \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_fiducial_csvs \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_fiducial_csvs \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_visualigned_imgs \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_qnii_nlaligned_coords

Created by Mukund on 2022-03-30

"""


from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import csv
import numpy as np
import json
import subprocess

nissl_id = int(sys.argv[1])
visualign_json =  sys.argv[2]
input_imgs = sys.argv[3]
chuck_space_coords_folder = sys.argv[4]
tfm_folder = sys.argv[5]
intrim_fiducial_folder = sys.argv[6]
op_imgs_folder = sys.argv[7]
op_qnii_coords_folder = sys.argv[8]

if (nissl_id<0):
    print("mapping all")
    # nissl_ids = [141, 143]
    # nissl_ids = [123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143]
    nissl_ids = [1, 3, 133, 135, 137, 139, 143]
    nissl_ids = [i for i in np.arange(1,214,2)]
    nissl_ids.remove(5)
    nissl_ids.remove(77)
    nissl_ids.remove(167)

    # nissl_ids.remove(181)
    # nissl_ids.remove(205)
    # nissl_ids.remove(223)
    # nissl_ids.remove(225)
    # nissl_ids.remove(227)
else:
    nissl_ids = [nissl_id]

print (nissl_id)


# read json and create csv fiducial files for cmapper

# read json
f = open(visualign_json)

# returns JSON object as
# a dictionary
qnii_data_unsorted = json.load(f)["slices"]

qnii_data = {}

for item in qnii_data_unsorted:
    try:
        nissl_id = int(item["filename"][4:7])
        markers = item["markers"]
        qnii_data[nissl_id] = markers
    except:
        print("VisuAlign transform missing for :", nissl_id)

    # img_dims[nissl_id] = {"height": item["height"], "width":item["width"]}
# print(qnii_data_unsorted)

for nissl_id in nissl_ids:
    # warped_coords
    nis_idx = str(nissl_id).zfill(3)

    # # convert fiducials to csv format
    # fids = np.asarray(qnii_data[nissl_id])
    # from_fids = fids[:,:][:, :2]
    # to_fids = fids[:,:][:,2:]
    # ids = np.array(range(len(fids)))
    # ids = ids.reshape((len(fids),1))
    # from_fids = np.concatenate((ids, from_fids), axis=1)
    # to_fids = np.concatenate((ids, to_fids), axis=1)

    # from_fiducials_file = intrim_fiducial_folder+"/"+str(nis_idx)+"_f.csv"
    # np.savetxt(from_fiducials_file, from_fids, delimiter=",")

    # to_fiducials_file = intrim_fiducial_folder+"/"+str(nis_idx)+"_t.csv"
    # np.savetxt(to_fiducials_file, to_fids, delimiter=",")

    # # tfm_folder = chuck_space_coords_folder
    # # intrim_pts_file = tfm_folder+"/"+nis_idx+"_intrim_pts.txt"
    # # with open(intrim_pts_file, 'w', newline='\n') as csvfile:
    # #     writer = csv.writer(csvfile)
    # #     for row in pts_tfmed_ss_nifti:
    # #         writer.writerow(row)

    intrim_pts_file = chuck_space_coords_folder+"/chuck_sp_img_coords_"+nis_idx+".csv" # pts to be transformed aka input pts

    # tfm_folder = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s6_register_to_allen/s6_fiducial_csvs"

    from_fiducials_file = tfm_folder+"/"+str(nissl_id)+"_f.csv"
    to_fiducials_file = tfm_folder+"/"+str(nissl_id)+"_t.csv"
    # nrrd_path = labelmap_folder+"/lmap_"+nis_idx+".nrrd"
    temp_output_csv_file = op_qnii_coords_folder+"/qnii_coords_"+nis_idx+".csv"
    # call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
    subprocess.run(["./build/cmapper2",
                    from_fiducials_file,
                    to_fiducials_file,
                    intrim_pts_file,
                    "None",
                    nis_idx,
                    temp_output_csv_file])
