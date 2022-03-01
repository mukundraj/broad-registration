""" 

Script to get region ids from puck positions. Brings togther various other
scripts to accomplish substeps.

Created by Mukund on 2022-02-16

Inputs: - path to folder with transform files ({slide_id}_f.csv {slide_id}_t.csv and {slide_id}_t3a.txt)
        - slide_idx (comes from nissl slide idx)
        - path to nrrd file with region ids
        - slice_idx (slice id in nrrd during sampling)
        - input_points csv file with puck positions to find segment id for
        - output.csv an output filename to store segment ids

Output: - Writes out csv with segment ids

Usage example: 

python src/python/scripts/get_region_ids.py \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/transforms \
79 \
/Users/mraj/Desktop/segmentations/79-83_143_labels.nrrd \
0 \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/testing/input_points.csv \
output/output.csv

References:

"""

import sys
import subprocess

# read in arguments
transforms_folder = sys.argv[1]
slide_idx = sys.argv[2]
nrrd_path = sys.argv[3]
slice_idx = sys.argv[4]
input_pts_file = sys.argv[5]
output_segids_file = sys.argv[6]

prealign_transform_path = transforms_folder+"/"+str(slide_idx)+"_t3a.txt"
intrim_t3a_pts_file = transforms_folder+"/"+str(slide_idx)+"_intrim_t3a.txt"

# call script to transform input points using prealignment transform
subprocess.run(["python", 
                "src/python/scripts/mapper_t3a.py", 
                prealign_transform_path, 
                input_pts_file, 
                intrim_t3a_pts_file])

from_fiducials_file = transforms_folder+"/"+str(slide_idx)+"_f.csv"
to_fiducials_file = transforms_folder+"/"+str(slide_idx)+"_t.csv"
# call cpp script to transform output of previous step with fiduals/TPS, sample from nrrd and write output
subprocess.run(["./build/cmapper", 
                from_fiducials_file,
                to_fiducials_file,
                intrim_t3a_pts_file,
                nrrd_path,
                str(slice_idx),
                output_segids_file])
