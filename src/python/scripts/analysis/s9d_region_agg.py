
"""
Generates data for regionally aggregated and normalized gene expression
visualization. Should be run after s9d_labelgen.py. Calls R script s9d_prep.R
internally for aggregating gene expression counts over CCF regions.

Usage:

python s9d_region_agg.py \
    inp: data root
    inp: path to processed region labels created by s9d_labelgen.py
    inp: path to full bead info including label info 
    inp: path to nissl images
    inp: path to wireframe images
    inp: path to chuck space coords
    i/o: path to interim folder for R script to write to 
    out: path to output folder

Usage example:

python src/python/scripts/analysis/s9d_region_agg.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/region_labels \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png \
    /data_v3_nissl_post_qc/s7_annotations/allen_labels_imgs/wireframe \
    /data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    1 3 \
    /data_v3_nissl_post_qc/s8_raw_data/integrated_mats \
    /data_v3_nissl_post_qc/s9_analysis/s9d/interim \
    /data_v3_nissl_post_qc/s9_analysis/s9d/gene_csvs_s9d

Created by Mukund on 2022-07-07

References:

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/gene_csvs_s9d gs://ml_portal2/test_data2/

"""

import anndata as ann
from scipy.sparse import csr_matrix, hstack
from produtils import dprint
import gc
import csv
import numpy as np
import os
import shutil
import json
import subprocess
import sys

data_root = sys.argv[1]
ip_folder_labels = data_root+sys.argv[2]
label_data_folder = data_root+sys.argv[3]
ip_folder_nissl = data_root+sys.argv[4]
ip_folder_atlas = data_root+sys.argv[5]
ip_folder_chuck_coords = data_root+sys.argv[6]
start_pid = sys.argv[7]
end_pid = sys.argv[8]
ip_folder_counts = data_root+sys.argv[9]
io_folder_interim = data_root+sys.argv[10]
op_folder = data_root+sys.argv[11]

# ip_folder_labels = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/region_labels'
# label_data_folder = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds'
# ip_folder_nissl = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png'
# ip_folder_atlas = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/allen_labels_imgs/wireframe'
# ip_folder_chuck_coords = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds'
# io_folder_interim = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/interim'
# op_folder = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/gene_csvs_s9d'

# file = "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/interim/output.h5ad"
# X = ann.read_h5ad(file)

# # create obs_names to idx map
# obs_names_list = list(X.obs_names)
# var_names_list = list(X.var_names)

# region_to_idx = {}
# for idx, obs in enumerate(obs_names_list):
#     region_to_idx[obs] = idx

# print(X.X)
# print(X.obs_names)
# print(X.var_names)
# csr = csr_matrix(X.X)
# dprint(csr.shape)

# print(csr.getrow(1))
subprocess.run(["Rscript", "src/python/scripts/analysis/s9d_prep.R", \
                ip_folder_labels, ip_folder_counts, io_folder_interim, \
                start_pid, end_pid])


pids = list(range(int(start_pid), int(end_pid)+1, 2))

for pid in pids:

    dprint(f'starting pid {pid}..................')
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    nis_id_str = str(apid).zfill(3)

    labels_csv_file = f'{label_data_folder}/allen_anno_data_{nis_id_str}.csv'
    # dprint(labels_csv_file)
    region_names = []
    out_tissue = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])
            out_tissue.append(row[8])

    puck_folder = f'{op_folder}/puck{pid}'
    if os.path.exists(puck_folder):
        shutil.rmtree(puck_folder)
    os.mkdir(puck_folder)


    # get chuck space img coords
    chuck_sp_img_coords_file = f'{ip_folder_chuck_coords}/chuck_sp_img_coords_{nis_id_str}.csv'
    chuck_sp_img_coords = []
    with open(chuck_sp_img_coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            chuck_sp_img_coords.append([int(row[0]), int(row[1])])

    # copy nissl and atlas images to puck directory
    from_nis_file = f'{ip_folder_nissl}/nis_{nis_id_str}.png'
    to_nis_file = f'{puck_folder}/nis_{nis_id_str}.png'
    from_atlas_file = f'{ip_folder_atlas}/chuck_sp_wireframe_{nis_id_str}.png'
    to_atlas_file = f'{puck_folder}/chuck_sp_wireframe_{nis_id_str}.png'
    # dprint(from_atlas_file)
    # dprint(to_atlas_file)
    subprocess.run(["cp", from_nis_file, to_nis_file])
    subprocess.run(["cp", from_atlas_file, to_atlas_file])


    # get labels for current puck's beads
    labels_csv_file = f'{ip_folder_labels}/agg_labels_{nis_id_str}.csv'
    dprint(labels_csv_file)
    region_ids = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_ids.append(int(row[0]))

    # writing coords tsv
    coords_csv_name = f'{puck_folder}/coords.csv'
    # dprint(np.shape(coords_dense_np))
    # dprint(coords_csv_name, pid, apid)
    # np.savetxt(coords_csv_name, np.array([xs,ys,zs]).T, fmt='%i', header="x,y,z", comments='', delimiter=",")
    in_tissue_region_ids = []
    in_tissue_inds = []
    with open(coords_csv_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=':')
        # writer.writerow(['x', 'y', 'z', 'rname'])
        writer.writerow(['x', 'y', 'rname'])
        for i, status in enumerate(out_tissue):
            if (status=='F'):
                # writer.writerow([xs[i], ys[i], zs[i], region_names[i]])
                writer.writerow([chuck_sp_img_coords[i][0], chuck_sp_img_coords[i][1], region_names[i]])
                in_tissue_inds.append(i)
                in_tissue_region_ids.append(region_ids[i])


    # get aggregated counts file
    aggr_counts_file = f'{io_folder_interim}/aggr_counts_{nis_id_str}.h5ad'
    aggr_counts = ann.read_h5ad(aggr_counts_file)

    # create obs_names to idx map
    obs_names_list = list(aggr_counts.obs_names)
    region_to_idx = {}
    for idx, obs in enumerate(obs_names_list):
        region_to_idx[int(obs)] = idx

    print(region_to_idx)

    nonzero_genes = []

    # iterate over genes
    genes = list(aggr_counts.var_names)
    for gene_idx, gene in enumerate(genes):
        # if gene=='Pcp4' or gene=='Gad2':
        #     dprint(f'Found {gene} at {gene_idx}')
        # else:
        #     continue

        if (gene_idx%500==0):
            collected = gc.collect()
            dprint('gene_idx', gene_idx, 'pid', pid, 'collected', collected)

        dprint(csr_matrix(aggr_counts.X).shape)

        spec_gene_regagg_cnts = csr_matrix(aggr_counts.X).getcol(gene_idx)
        spec_gene_regagg_cnts_dense = np.squeeze(np.array(spec_gene_regagg_cnts.todense())).astype(int)
        dprint(len(in_tissue_inds))
        dprint(spec_gene_regagg_cnts_dense.shape)

        # read in bead counts and normalize
        ip_bead_counts_file = f'{io_folder_interim}/aggr_num_beads_{nis_id_str}.csv'
        bead_counts = np.genfromtxt(ip_bead_counts_file, delimiter=',', skip_header=True, dtype=np.int32)
        bead_counts = bead_counts[:,1]

        # dprint(spec_gene_regagg_cnts_dense)
        dprint(np.shape(spec_gene_regagg_cnts_dense))
        spec_gene_regagg_cnts_dense = np.squeeze(spec_gene_regagg_cnts_dense/bead_counts[None,:])
        dprint(np.shape(np.squeeze(spec_gene_regagg_cnts_dense)))

        # dprint(bead_counts)
        regional_aggr_gene_cnts = np.zeros(len(in_tissue_region_ids))

        dprint(np.shape(spec_gene_regagg_cnts_dense))
        # iterate over beads in puck and update aggregated counts
        for idx in range(len(regional_aggr_gene_cnts)):
            label = in_tissue_region_ids[idx]
            reg_idx = region_to_idx[label]
            regional_aggr_gene_cnts[idx] = spec_gene_regagg_cnts_dense[reg_idx]

        # check if gene expressed at all in any region in puck
        if (np.sum(regional_aggr_gene_cnts)>0):
            nonzero_genes.append(gene)

            # write out counts and metadata
            gene_metadata={"maxCount":np.max(regional_aggr_gene_cnts)}
            gene_cnts=regional_aggr_gene_cnts
            gene_csv_name = f'{puck_folder}/gene_{gene}.csv'
            np.savetxt(gene_csv_name, gene_cnts, fmt='%.2e', header="count", comments='',delimiter=',')

            metadata_json_file = f'{puck_folder}/metadata_gene_{gene}.json'
            with open(metadata_json_file, 'w') as outfile:
                tmp_dict = {'maxCount':str(gene_metadata['maxCount'])}
                json.dump(tmp_dict, outfile, separators=(',', ':'))

    geneOptions_json_file = f'{puck_folder}/geneOptions.json'
    gene_options_dict = {'geneOptions':nonzero_genes}
    with open(geneOptions_json_file, 'w') as outfile:
        json.dump(gene_options_dict, outfile, separators=(',', ':'))

