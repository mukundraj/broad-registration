"""Generate new labelmap csvs based on labels derived from Partha contours

Usage:

python s2b_gen_labelmap_csvs.py \
    inp: data_root \
    inp: path to bead coords in chuck space
    inp: path to anndata files with new Partha labels
    inp: name to acro map
    out: path to save new labelmap csvs:w

Example:

python src/python/scripts/v3/s2b_gen_labelmap_csvs.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /v3/s1/CSHL_CTMapping \
    /v3/s1/allen_name_to_acro_map.csv \
    /v3/s2/bead_ccf_labels_allbds \

Created by Mukund on 2022-12-22

"""

import sys
import csv
from produtils import dprint
import anndata as ann


data_root = sys.argv[1]
ip_folder_chuck_coords = f'{data_root}{sys.argv[2]}'
anndata_path = f'{data_root}{sys.argv[3]}'
name_to_acro_map_file = f'{data_root}{sys.argv[4]}'
new_labelmap_csvs_path = f'{data_root}{sys.argv[5]}'

# read in allen_name_to_acro_map
name_to_id_map = {}
name_to_acro_map = {}
with open(name_to_acro_map_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        name_to_acro_map[row[2]] = row[1]
        name_to_id_map[row[2]] = row[0]

name_to_acro_map['NA'] = 'NA'
name_to_id_map['NA'] = '0'
name_to_id_map['OUT'] = '0'


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
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    dprint('starting pid', pid)
    # get chuck space img coords
    chuck_sp_img_coords_file = f'{ip_folder_chuck_coords}/chuck_sp_img_coords_{nis_id_str}.csv'
    chuck_sp_img_coords = []
    with open(chuck_sp_img_coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            chuck_sp_img_coords.append([int(row[0]), int(row[1])])

    # create new labelmap array with default OUT value
    nbeads = len(chuck_sp_img_coords)
    labels = ['OUT']*nbeads
    outTissueStatus = ['T']*nbeads
    acronyms = ['NA']*nbeads

    # read anndata file with partha labels
    anndata_file = f'{anndata_path}/dd{apid}_CTMapping.h5ad'
    counts = ann.read_h5ad(anndata_file)
    position_inds = counts.var.rowName.values.astype(int)
    cshl_labels = counts.var.CSHL_CCFname.values


    # iterate over beads and update labels for IN beads
    for idx, bidx in enumerate(position_inds):
        zbidx = bidx # zero based bead index
        labels[zbidx] = cshl_labels[idx]
        outTissueStatus[zbidx] = 'F'
        # dprint(zbidx, labels[zbidx], cshl_labels[idx])
        acronyms[zbidx] = name_to_acro_map[labels[zbidx]]

    # save new labelmap csv
    new_labelmap_csv_file = f'{new_labelmap_csvs_path}/allen_anno_data_{nis_id_str}.csv'
    with open(new_labelmap_csv_file, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        for idx in range(nbeads):
            writer.writerow(['','','',acronyms[idx], labels[idx],'','','',outTissueStatus[idx], name_to_id_map[labels[idx]]])

