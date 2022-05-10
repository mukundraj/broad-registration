""" Script to process annodata files generated by R script (which generates the
annodata from Jonah's qs object). Output expected to be zarr file to be used
downstream for interactive visualization.

Usage:

python s8_process_integ_mat.py
    inp: input folder
    out: intermediate zarr dir
    out: final zarr dir

Usage example:

python src/python/scripts/v2/s8_process_integ_mat.py \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/gene_jsons

Created by Mukund on 2022-05-04

References:
- on sparse and zarr https://github.com/zarr-developers/zarr-python/issues/424

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

in_folder = sys.argv[1]
op_folder = sys.argv[2]

genes_list = ['Pcp4', 'Calb1', 'Gng13', 'Gabra6']

for pid in range(1,4):

    ip_coords_file  = f'{in_folder}/ad_coords_{str(pid)}.h5ad'
    ip_counts_file  = f'{in_folder}/ad_counts_{str(pid)}.h5ad'

    counts = ann.read_h5ad(ip_counts_file)
    coords = ann.read_h5ad(ip_coords_file)

    counts_X = csr_matrix(counts.X).transpose()
    coords_X = csr_matrix(coords.X)
    coords_dense_np = np.array(coords_X.todense())
    xs = coords_dense_np[:, 0].astype(int).tolist()
    ys = coords_dense_np[:, 1].astype(int).tolist()
    zs = coords_dense_np[:, 2].astype(int).tolist()
    data = {'x': xs,
            'y': ys,
            'z': zs}
    json_string = json.dumps(data)

    puck_folder = f'{op_folder}/puck{pid}'
    if os.path.exists(puck_folder):
        shutil.rmtree(puck_folder)
    os.mkdir(puck_folder)

    # writing coords tsv
    coords_csv_name = f'{puck_folder}/coords.csv'
    dprint(np.shape(coords_dense_np))
    dprint(coords_csv_name)
    np.savetxt(coords_csv_name, np.array([xs,ys,zs]).T, fmt='%i', header="x,y,z", comments='', delimiter=",")

    json_file = f'{puck_folder}/coords.json'
    # Directly from dictionary
    with open(json_file, 'w') as outfile:
        json.dump(data, outfile, separators=(',', ':'))

    genes = list(counts.obs_names)

    gene_cnts = {}
    for gene in genes_list:
        gene_idx = genes.index(gene)
        specific_gene_cnts = counts_X.getcol(gene_idx)
        spec_gene_cnts_dense = np.squeeze(np.array(specific_gene_cnts.todense())).astype(int)
        gene_cnts[gene]=spec_gene_cnts_dense

    for key in gene_cnts:
        gene_csv_name = f'{puck_folder}/gene_{key}.csv'
        np.savetxt(gene_csv_name, gene_cnts[key], fmt='%i', header="count", comments='',delimiter=',')
        json_file = f'{puck_folder}/gene_{key}.json'
        # Directly from dictionary
        with open(json_file, 'w') as outfile:
            tmp_dict = {key:json.dumps(gene_cnts[key].tolist())}
            json.dump(tmp_dict, outfile, separators=(',', ':'))

    dprint(f'puck {pid} done')

