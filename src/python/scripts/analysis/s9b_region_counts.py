"""
Script for creating matrices containing gene counts per region (unnormalized and normalized).

Usage:

python s9_region_counts.py \
    inp: path to bead_ccf_labels_allbds with bead metadata in csv format
    inp: path to bead_ccf_coords_allbds with bead coords in chuck space in csv format
    inp: input to integrated_mats folder with processed gene counts in annodata h5ad format
    inp: file (json) with region names in alphabetical order with associated index
    inp: path to nissl images folder
    inp: path to atlas wireframe images folder
    out: output dir to write gene_csvs

Usage example:

python src/python/scripts/analysis/s9b_region_counts.py \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/ccf_regions.json \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/allen_labels_imgs/wireframe \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/gene_csvs \

Created by Mukund on 2022-05-19

References - https://stackoverflow.com/questions/9001509/how-can-i-sort-a-dictionary-by-key

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/gene_csvs gs://ml_portal2/test_data2/

"""
from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))

import src.python.utils.io as io

import sys
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix, hstack
import zarr
import json
from produtils import dprint
import os
import shutil
import csv
import subprocess

ip_folder_labels = sys.argv[1]
ip_folder_chuck_coords = sys.argv[2]
ip_folder_counts = sys.argv[3]
ip_file_ccf_regions = sys.argv[4]
ip_folder_nissl = sys.argv[5]
ip_folder_atlas = sys.argv[6]
op_folder_gene_csvs = sys.argv[7]

genes_list = ['Gad1', 'Gad2', 'Slc17a7']

region_names_dict = {}

with open(ip_file_ccf_regions) as json_file:
    region_names_dict = json.load(json_file)

# dprint(region_names_dict)
nRegions = len(region_names_dict.keys())
dprint("nGenes : ", nRegions)

# get number of regions
for pid in range(1, 208, 2):

    dprint(f'starting pid {pid}..................')
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # ip_coords_file  = f'{ip_folder_counts}/ad_coords_{str(apid)}.h5ad'
    ip_counts_file  = f'{ip_folder_counts}/ad_counts_{str(apid)}.h5ad'

    # get region names/labels
    nis_id_str = str(apid).zfill(3)
    labels_csv_file = f'{ip_folder_labels}/allen_anno_data_{nis_id_str}.csv'
    dprint(labels_csv_file)
    region_names = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])


    # get indices of beads on puck inside tissue region, write out coords to destination folder
    puck_folder = io.recreate_puck_folder(op_folder_gene_csvs, apid)
    in_tissue_inds = io.read_mask_write_beads(apid, ip_folder_labels, ip_folder_chuck_coords, ip_folder_counts, op_folder_gene_csvs)
    dprint(len(in_tissue_inds))

    # copy nissl and atlas images to puck directory
    from_nis_file = f'{ip_folder_nissl}/nis_{nis_id_str}.png'
    to_nis_file = f'{puck_folder}/nis_{nis_id_str}.png'
    # from_atlas_file = f'{ip_folder_atlas}/chuck_sp_labelmap_{nis_id_str}.png'
    from_atlas_file = f'{ip_folder_atlas}/chuck_sp_wireframe_{nis_id_str}.png'
    to_atlas_file = f'{puck_folder}/chuck_sp_wireframe_{nis_id_str}.png'
    dprint(from_atlas_file)
    dprint(to_atlas_file)
    subprocess.run(["cp", from_nis_file, to_nis_file])
    subprocess.run(["cp", from_atlas_file, to_atlas_file])


    counts = ann.read_h5ad(ip_counts_file)
    counts_X = csr_matrix(counts.X).transpose()
    genes = list(counts.obs_names)

    # for each gene
    for gene in genes_list:

        gene_idx = genes.index(gene)
        specific_gene_cnts = counts_X.getcol(gene_idx)
        spec_gene_cnts_dense = np.squeeze(np.array(specific_gene_cnts.todense())).astype(int)

        region_beads = np.zeros(nRegions)
        region_counts = np.zeros(nRegions)
        # for each bead in tissue, increase count toward respective region
        for idx in in_tissue_inds:
            region_name = region_names[idx]
            region_idx = region_names_dict[region_name]
            region_beads[region_idx] += 1
            region_counts[region_idx] += spec_gene_cnts_dense[idx]

        bead_counts = np.zeros(len(in_tissue_inds))
        bead_normed_counts = np.zeros(len(in_tissue_inds))
        # for each in_tissue bead assign gene's count + normed count in its region
        for i, idx in enumerate(in_tissue_inds):
            region_name = region_names[idx]
            region_idx = region_names_dict[region_name]

            bead_counts[i] = region_counts[region_idx]
            bead_normed_counts[i] = float(region_counts[region_idx])/region_beads[region_idx]

        # write out gene csv file with total regional count for each bead
        gene_csv_name = f'{puck_folder}/rc_{gene}.csv'
        np.savetxt(gene_csv_name, bead_counts, fmt='%i', header="count", comments='',delimiter=',')
        # json_file = f'{puck_folder}/gene_{gene}.json'


        # write out gene csv file with normalized regional count for each bead
        gene_csv_name = f'{puck_folder}/rnc_{gene}.csv'
        np.savetxt(gene_csv_name, bead_normed_counts, fmt='%f', header="count", comments='',delimiter=',')
        # json_file = f'{puck_folder}/gene_{gene}.json'

        max_count = np.max(bead_counts)
        max_normed_count = np.max(bead_normed_counts)
        metadata_json_file = f'{puck_folder}/metadata_gene_{gene}.json'
        with open(metadata_json_file, 'w') as outfile:
            tmp_dict = {'maxCount':str(max_count), 'maxNormedCount':str(round(max_normed_count,3))}
            json.dump(tmp_dict, outfile, separators=(',', ':'))


        dprint(bead_counts)
        dprint(bead_normed_counts)
        # gene_cnts = {}
        # gene_metadata = {}
        # for gene in genes_list:
        #     gene_idx = genes.index(gene)
        #     specific_gene_cnts = counts_X.getcol(gene_idx)
        #     spec_gene_cnts_dense = np.squeeze(np.array(specific_gene_cnts.todense())).astype(int)
        #     spec_gene_cnts_dense = spec_gene_cnts_dense[in_tissue_inds]
        #     dprint(np.max(spec_gene_cnts_dense))
        #     gene_metadata[gene]={"maxCount":np.max(spec_gene_cnts_dense)}
        #     gene_cnts[gene]=spec_gene_cnts_dense

        # for key in gene_cnts:
        #     gene_csv_name = f'{puck_folder}/gene_{key}.csv'
        #     np.savetxt(gene_csv_name, gene_cnts[key], fmt='%i', header="count", comments='',delimiter=',')
        #     json_file = f'{puck_folder}/gene_{key}.json'
        #     # Directly from dictionary
        #     with open(json_file, 'w') as outfile:
        #         tmp_dict = {key:json.dumps(gene_cnts[key].tolist())}
        #         json.dump(tmp_dict, outfile, separators=(',', ':'))

        #     metadata_json_file = f'{puck_folder}/metadata_gene_{key}.json'
        #     with open(metadata_json_file, 'w') as outfile:
        #         dprint(gene_metadata[key])
        #         dprint(key)
        #         tmp_dict = {'maxCount':str(gene_metadata[key]['maxCount'])}
        #         json.dump(tmp_dict, outfile, separators=(',', ':'))



