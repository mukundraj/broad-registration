"""
Script for generating identifing and storing bead indices from specific regions

Usage:

python s9c_qc_prep.py \
    inp: data root path
    inp: path to bead_ccf_labels_allbds with bead metadata in csv format
    inp: path to bead_ccf_coords_allbds with bead coords in chuck space in csv format
    inp: input to integrated_mats folder with processed gene counts in annodata h5ad format
    out: output dir to write qc_mats

Usage example:

python src/python/scripts/analysis/s9c_qc_prep.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats \
    /data_v3_nissl_post_qc/s9_analysis/qc_mats

Created by Mukund on 2022-06-09
"""
import sys
from produtils import dprint
import csv
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix, hstack
import json

from allensdk.core.reference_space_cache import ReferenceSpaceCache

reference_space_key = 'annotation/ccf_2017'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1) 

ip_data_root = sys.argv[1]
ip_folder_labels = ip_data_root+sys.argv[2]
ip_folder_chuck_coords = ip_data_root+sys.argv[3]
ip_folder_counts = sys.argv[4]
op_folder = ip_data_root+sys.argv[5]

# regional_inds = {
#     "Primary somatosensory area, trunk, layer 1":[],
#     "Primary somatosensory area, trunk, layer 2/3":[],
#     "Primary somatosensory area, trunk, layer 4":[],
#     "Primary somatosensory area, trunk, layer 5":[],
#     "Primary somatosensory area, trunk, layer 6a":[],
#     "Primary somatosensory area, trunk, layer 6b":[],
#     "Hippocampal formation":[]
# }

region_name_to_id = {
    'Primary somatosensory area, trunk, layer 1':-1,
    'Primary somatosensory area, trunk, layer 2/3':-1,
    'Primary somatosensory area, trunk, layer 4':-1,
    'Primary somatosensory area, trunk, layer 5':-1,
    'Primary somatosensory area, trunk, layer 6a':-1,
    'Primary somatosensory area, trunk, layer 6b':-1,
    'Hippocampal formation':-1
}
for rname in region_name_to_id:
    region_name_to_id[rname] = tree.get_structures_by_name([rname])[0]['id']
dprint(region_name_to_id)

genes_list = ['Dcn', 'Ndnf']

pucks = list(range(1, 208, 2))
data_dicts = []
dprint(pucks)
for i in range(len(genes_list)):
    data = {"x": list(range(1, len(pucks)+1)), "y": [1,2], "z":[[],[]]}
    data_dicts.append(data)
dprint(data_dicts)
# exit(0)
# get number of regions
for pid in pucks:

    dprint("starting pid ", pid)
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # ip_coords_file  = f'{ip_folder_counts}/ad_coords_{str(apid)}.h5ad'
    ip_counts_file  = f'{ip_folder_counts}/ad_counts_{str(apid)}.h5ad'

    # get region names/labels
    nis_id_str = str(apid).zfill(3)
    labels_csv_file = f'{ip_folder_labels}/allen_anno_data_{nis_id_str}.csv'
    region_names = []
    out_tissue = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])
            out_tissue.append(row[8])

    counts = ann.read_h5ad(ip_counts_file)
    counts_X = csr_matrix(counts.X).transpose()
    genes = list(counts.obs_names)
    gene_count_total = counts_X.sum(axis=1)
    # dprint(counts_X.get_shape())
    dprint('gene_count_total shape:', np.shape(gene_count_total))

    nis_id_str = str(apid).zfill(3)
    # get chuck space img coords
    chuck_sp_img_coords_file = f'{ip_folder_chuck_coords}/chuck_sp_img_coords_{nis_id_str}.csv'
    dprint('chuck_sp_img_coords file ', chuck_sp_img_coords_file)
    chuck_sp_img_coords = []
    with open(chuck_sp_img_coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            chuck_sp_img_coords.append([int(row[0]), int(row[1])])
    dprint(len(chuck_sp_img_coords))
    max_len = len(chuck_sp_img_coords)

    regional_inds = {
        "Primary somatosensory area, trunk, layer 1":[],
        "Primary somatosensory area, trunk, layer 2/3":[],
        "Primary somatosensory area, trunk, layer 4":[],
        "Primary somatosensory area, trunk, layer 5":[],
        "Primary somatosensory area, trunk, layer 6a":[],
        "Primary somatosensory area, trunk, layer 6b":[],
        "Hippocampal formation":[]
    }

    for idx in range(len(chuck_sp_img_coords)):
        leaf_region_name = region_names[idx]
        assert(idx<max_len)
        if (out_tissue[idx]=='F'):
            # if (idx>170950):
            #     dprint(idx)
            leaf_id = tree.get_structures_by_name([leaf_region_name])[0]['id']
            # dprint(leaf_region_name, leaf_id)
            # if tree.structure_descends_from(leaf_id, 1006):
            if leaf_id==1006:
                # Primary somatosensory area, trunk, layer 1
                regional_inds['Primary somatosensory area, trunk, layer 1'].append(idx)
            # elif tree.structure_descends_from(leaf_id, 670):
            elif leaf_id==670:
                # Primary somatosensory area, trunk, layer 2/3
                regional_inds['Primary somatosensory area, trunk, layer 2/3'].append(idx)
            # elif tree.structure_descends_from(leaf_id, 1086):
            elif leaf_id==1086:
                # Primary somatosensory area, trunk, layer 4
                regional_inds['Primary somatosensory area, trunk, layer 4'].append(idx)
            # elif tree.structure_descends_from(leaf_id, 1111):
            elif leaf_id==1111:
                # Primary somatosensory area, trunk, layer 5
                regional_inds['Primary somatosensory area, trunk, layer 5'].append(idx)
            # elif tree.structure_descends_from(leaf_id, 9):
            elif leaf_id==9:
                # Primary somatosensory area, trunk, layer 6a
                regional_inds['Primary somatosensory area, trunk, layer 6a'].append(idx)
            # elif tree.structure_descends_from(leaf_id, 461):
            elif leaf_id==461:
                # Primary somatosensory area, trunk, layer 6b
                regional_inds['Primary somatosensory area, trunk, layer 6b'].append(idx)
            elif tree.structure_descends_from(leaf_id, 1089):
                # Hippocampal formation
                regional_inds['Hippocampal formation'].append(idx)

            # if tree.structure_descends_from(str_id, cerebral_cortex_id):
            #     gene_count_cortex += spec_gene_cnts_dense[idx]
            #     total_count_cortex += gene_count_total[idx][0,0]
            # if tree.structure_descends_from(str_id, hippo_id):
            #     gene_count_hippo += spec_gene_cnts_dense[idx]
            #     total_count_hippo += gene_count_total[idx][0,0]
    # dprint(regional_inds)
    with open(f'{op_folder}/regional_inds_{pid}.json', "w") as outfile:
        json.dump(regional_inds, outfile)
