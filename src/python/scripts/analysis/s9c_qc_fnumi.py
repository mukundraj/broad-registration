"""
Script for generating data for QC heatmap with fractional nUMIs (gene X over all nUMIs) 

Usage:

python s9c_qc_fnumi.py \
    inp: data root path
    inp: path to bead_ccf_labels_allbds with bead metadata in csv format
    inp: path to bead_ccf_coords_allbds with bead coords in chuck space in csv format
    inp: input to integrated_mats folder with processed gene counts in annodata h5ad format
    out: output dir to write qc_mats

Usage example:

python src/python/scripts/analysis/s9c_qc_fnumi.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats \
    /data_v3_nissl_post_qc/s9_analysis/qc_mats

Created by Mukund on 2022-06-08

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

    # for each gene
    for gid, gene in enumerate(genes_list):

        gene_idx = genes.index(gene)
        specific_gene_cnts = counts_X.getcol(gene_idx)
        spec_gene_cnts_dense = np.squeeze(np.array(specific_gene_cnts.todense())).astype(int)
        # dprint(np.shape(spec_gene_cnts_dense))

        nis_id_str = str(apid).zfill(3)
        # get chuck space img coords
        chuck_sp_img_coords_file = f'{ip_folder_chuck_coords}/chuck_sp_img_coords_{nis_id_str}.csv'
        chuck_sp_img_coords = []
        with open(chuck_sp_img_coords_file, newline='\n') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                chuck_sp_img_coords.append([int(row[0]), int(row[1])])

        # dprint(len(chuck_sp_img_coords))

        total_count_cortex = 0
        total_count_hippo = 0
        gene_count_cortex = 0
        gene_count_hippo = 0
        cerebral_cortex_id = tree.get_structures_by_name(['Cerebral cortex'])[0]['id']
        hippo_id = tree.get_structures_by_name(['Hippocampal formation'])[0]['id']
        dprint('cerid', cerebral_cortex_id, 'hippoid:', hippo_id)
        for idx in range(len(chuck_sp_img_coords)):
            region_name = region_names[idx]
            if (out_tissue[idx]=='F'):
                str_id = tree.get_structures_by_name([region_name])[0]['id']
                if tree.structure_descends_from(str_id, cerebral_cortex_id):
                    gene_count_cortex += spec_gene_cnts_dense[idx]
                    total_count_cortex += gene_count_total[idx][0,0]
                if tree.structure_descends_from(str_id, hippo_id):
                    gene_count_hippo += spec_gene_cnts_dense[idx]
                    total_count_hippo += gene_count_total[idx][0,0]
                # dprint(idx,str_id)

            # region_idx = region_names_dict[region_name]
            # region_beads[region_idx] += 1
            # region_counts[region_idx] += spec_gene_cnts_dense[idx]
        dprint(gene_count_cortex, gene_count_hippo, total_count_cortex, total_count_hippo)
        if(gene_count_cortex>0):
            data_dicts[gid]['z'][0].append(float(gene_count_cortex)/total_count_cortex)
        else:
            data_dicts[gid]['z'][0].append(0)

        if (gene_count_hippo>0):
            data_dicts[gid]['z'][1].append(float(gene_count_hippo)/total_count_hippo)
        else:
            data_dicts[gid]['z'][1].append(0)
        dprint(data_dicts)
        # thoughts - div by zero error, is total count cortex/hippo out of inner loop possible?


analysis_metadata = []
for idx, meta_data_item in enumerate(data_dicts):
    analysis_metadata.append({"name": genes_list[idx], "filename":f'{genes_list[idx]}.json'})
    with open(f'{op_folder}/{analysis_metadata[idx]["filename"]}', "w") as outfile:
        json.dump(data_dicts[idx], outfile)

meta_data={"analysis_metadata":analysis_metadata}
with open(op_folder+"/metadata.json", "w") as outfile:
    json.dump(meta_data, outfile)
