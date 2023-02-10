"""Generate coords csv file with x,y in Chuck space and region id based on
Partha's contours. CAVEAT - Currently this script generates coords.csv files
only for GeneExp tab. Coords for CellSpatial tab are generaged in
s1c_beadxcell_zarr.py

Usage:

python s2c_gen_coords.py \
    inp: data_root \
    inp: path to bead_ccf_coords_allbds with bead coords in chuck space in csv format \
    inp: path to bead_ccf_coords_allbds with partha labels
    out: file path to save coords in chuck space with partha region id


Example:

python src/python/scripts/v3/s2c_gen_coords.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /v3/s2/bead_ccf_labels_allbds \
    /v3/s2/coords \



Created by Mukund on 2022-12-22
"""

from produtils import dprint
import sys
import csv


data_root = sys.argv[1]
ip_folder_chuck_coords = f'{data_root}{sys.argv[2]}'
ip_folder_partha_labels = f'{data_root}{sys.argv[3]}'

dprint(data_root)


start_pid = 1
end_pid = 207
pids = list(range(int(start_pid), int(end_pid)+1, 2))
if 5 in pids:
    pids.remove(5)
if 77 in pids:
    pids.remove(77)
if 167 in pids:
    pids.remove(167)
# iterate over pids
for pids_idx, pid in enumerate(pids):
    assert(pid!=5 and pid!=77 and pid!=167)
    nis_id_str = str(pid).zfill(3)
    # apid = pid
    # if (pid==5 or pid==77 or pid==167):
    #     apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # read in bead coords in chuck space
    chuck_sp_img_coords_file = f'{ip_folder_chuck_coords}/chuck_sp_img_coords_{nis_id_str}.csv'
    chuck_sp_img_coords = []
    with open(chuck_sp_img_coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            chuck_sp_img_coords.append([int(row[0]), int(row[1])])


    # read in partha labels
    partha_labels_file = f'{ip_folder_partha_labels}/allen_anno_data_{nis_id_str}.csv'
    partha_labels = []
    out_tissue = []
    with open(partha_labels_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            partha_labels.append(row[4])
            out_tissue.append(row[8])



    coords_csv_name = f'{data_root}/v3/s2/coords/coords_{pid}.csv'
    # generate celltype coords file
    with open(coords_csv_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=':')
        # writer.writerow(['x', 'y', 'z', 'rname'])
        writer.writerow(['x', 'y', 'rname'])
        for i, status in enumerate(out_tissue):
            if (status=='F'):
                # writer.writerow([xs[i], ys[i], zs[i], region_names[i]])
                writer.writerow([chuck_sp_img_coords[i][0], chuck_sp_img_coords[i][1], partha_labels[i]])
