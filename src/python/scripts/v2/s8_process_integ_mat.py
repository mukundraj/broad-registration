""" Script to process annodata files generated by R script (which generates the
annodata from Jonah's qs object). Output expected to be zarr file to be used
downstream for interactive visualization.

Usage:

python s8_process_integ_mat.py
    inp: path to bead_ccf_labels_allbds with bead metadata in csv format
    inp: path to bead_ccf_coords_allbds with bead coords in chuck space in csv format
    inp: input to integrated_mats folder with processed gene counts in annodata h5ad format
    inp: path to nissl images folder
    inp: path to atlas wireframe images folder
    out: output dir to write gene jsons to

Usage example:

python src/python/scripts/v2/s8_process_integ_mat.py \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/allen_labels_imgs/wireframe \
    /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/gene_jsons

Created by Mukund on 2022-05-04

References:
- on sparse and zarr https://github.com/zarr-developers/zarr-python/issues/424

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/gene_jsons gs://ml_portal2/test_data2/

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/gene_jsons/puck89 gs://ml_portal2/test_data2/gene_jsons/

"""

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

label_data_folder = sys.argv[1]
ip_folder_chuck_coords = sys.argv[2]
in_folder = sys.argv[3]
ip_folder_nissl = sys.argv[4]
ip_folder_atlas = sys.argv[5]
op_folder = sys.argv[6]

genes_list = ['Pcp4', 'Calb1', 'Gng13', 'Gabra6',
              'Mbp', 'Plp1', 'Mag',
              'Myoc', 'Agt', 'Gfap', 'Slc1a3', 'Aqp4',
              'Dcn', 'Flt1',
              'Rarres2', 'Foxj1']
genes_list.extend(['Xpo7', 'Cul1', 'Herc1', 'Rb1cc1', 'Setd1a','Trio', 
                   'Cacna1g', 'Sp4', 'Gria3', 'Grin2a',
                   'Slc17a7', 'Dsp','Gad2'])

# genes_list = ['slc17a7', 'Dsp', 'Gad2']
# genes_list = ['Pcp4']

# for pid in range(1,42,2):
for pid in range(1,208,2):

    dprint(f'starting pid {pid}..................')
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # ip_coords_file  = f'{in_folder}/ad_coords_{str(apid)}.h5ad'
    ip_counts_file  = f'{in_folder}/ad_counts_{str(apid)}.h5ad'

    counts = ann.read_h5ad(ip_counts_file)
    # coords = ann.read_h5ad(ip_coords_file)

    counts_X = csr_matrix(counts.X).transpose()
    # coords_X = csr_matrix(coords.X)
    # coords_dense_np = np.array(coords_X.todense())
    # xs = coords_dense_np[:, 0].astype(int).tolist()
    # ys = coords_dense_np[:, 1].astype(int).tolist()
    # zs = coords_dense_np[:, 2].astype(int).tolist()
    # data = {'x': xs,
    #         'y': ys,
    #         'z': zs}
    # json_string = json.dumps(data)

    # reading csv for label data

    nis_id_str = str(apid).zfill(3)
    labels_csv_file = f'{label_data_folder}/allen_anno_data_{nis_id_str}.csv'
    dprint(labels_csv_file)
    region_names = []
    out_tissue = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])
            out_tissue.append(row[8])

    dprint(len(region_names), len(out_tissue))
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

    # writing coords tsv
    coords_csv_name = f'{puck_folder}/coords.csv'
    # dprint(np.shape(coords_dense_np))
    dprint(coords_csv_name, pid, apid)
    # np.savetxt(coords_csv_name, np.array([xs,ys,zs]).T, fmt='%i', header="x,y,z", comments='', delimiter=",")
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

    # json_file = f'{puck_folder}/coords.json'
    # Directly from dictionary
    # with open(json_file, 'w') as outfile:
    #     json.dump(data, outfile, separators=(',', ':'))

    genes = list(counts.obs_names)

    gene_cnts = {}
    gene_metadata = {}
    for gene in genes_list:
        gene_idx = genes.index(gene)
        specific_gene_cnts = counts_X.getcol(gene_idx)
        spec_gene_cnts_dense = np.squeeze(np.array(specific_gene_cnts.todense())).astype(int)
        spec_gene_cnts_dense = spec_gene_cnts_dense[in_tissue_inds]
        dprint(np.max(spec_gene_cnts_dense))
        gene_metadata[gene]={"maxCount":np.max(spec_gene_cnts_dense)}
        gene_cnts[gene]=spec_gene_cnts_dense

    for key in gene_cnts:
        gene_csv_name = f'{puck_folder}/gene_{key}.csv'
        np.savetxt(gene_csv_name, gene_cnts[key], fmt='%i', header="count", comments='',delimiter=',')
        # json_file = f'{puck_folder}/gene_{key}.json'
        # # Directly from dictionary
        # with open(json_file, 'w') as outfile:
        #     tmp_dict = {key:json.dumps(gene_cnts[key].tolist())}
        #     json.dump(tmp_dict, outfile, separators=(',', ':'))

        metadata_json_file = f'{puck_folder}/metadata_gene_{key}.json'
        with open(metadata_json_file, 'w') as outfile:
            dprint(gene_metadata[key])
            dprint(key)
            tmp_dict = {'maxCount':str(gene_metadata[key]['maxCount'])}
            json.dump(tmp_dict, outfile, separators=(',', ':'))

    dprint(f'puck {apid} done')

