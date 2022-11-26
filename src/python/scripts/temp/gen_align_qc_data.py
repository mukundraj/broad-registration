""" Quickly put together script for generating data for Evan's alignment QC tests

python gen_align_qc_data.py
   i/p: data root folder
   inp: integrated mats
   inp: label data folder
   out: output folder

python src/python/scripts/temp/temp.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s8_raw_data/integrated_mats \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /scratch/gene_exp \


Created by Mukund on 2022-11-23

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
import time
from multiprocessing import Pool
import gc

data_root = sys.argv[1]
in_folder = data_root+sys.argv[2]
label_data_folder = data_root+sys.argv[3]
op_folder = data_root+sys.argv[4]


# for each puck read anndata

# get data for chosen genes 


# write csv with anndata
genes_list = ['Gad1', 'Gad2', 'Slc17a6', 'Slc17a7', 'Slc17a8', 'Rln3', 'Slc6a3', 'Slc6a4', 'Slc6a2', 'Dbh', 'Tph2', 'Th', 'Slc18a3', 'Chat', 'Dsp', 'Gm5741', 'Sln']

# genes_list = ['slc17a7', 'Dsp', 'Gad2']
# genes_list = ['Pcp4']

# for pid in range(1,42,2):
def process_pid(pid):

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
    # dprint(labels_csv_file)
    region_names = []
    out_tissue = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])
            out_tissue.append(row[8])

    # dprint(len(region_names), len(out_tissue))
    puck_folder = f'{op_folder}/puck{pid}'
    if os.path.exists(puck_folder):
        shutil.rmtree(puck_folder)
    os.mkdir(puck_folder)

    genes = list(counts.obs_names)

    for gene_idx, gene in enumerate(genes):
        # if gene=='Pcp4':
        if gene in genes_list:
            if (gene_idx%500==0):
                collected = gc.collect()
                dprint('gene_idx', gene_idx, 'pid', pid, 'collected', collected)
            specific_gene_cnts = counts_X.getcol(gene_idx)
            spec_gene_cnts_dense = np.squeeze(np.array(specific_gene_cnts.todense())).astype(int)
            # spec_gene_cnts_dense = spec_gene_cnts_dense[in_tissue_inds]
            # dprint(np.max(spec_gene_cnts_dense))
            # gene_metadata={"maxCount":np.max(spec_gene_cnts_dense)}
            gene_cnts=spec_gene_cnts_dense
            gene_csv_name = f'{puck_folder}/gene_{gene}.csv'
            np.savetxt(gene_csv_name, gene_cnts, fmt='%i', header="count", comments='',delimiter=',')

            # metadata_json_file = f'{puck_folder}/metadata_gene_{gene}.json'
            # with open(metadata_json_file, 'w') as outfile:
            #     tmp_dict = {'maxCount':str(gene_metadata['maxCount'])}
            #     json.dump(tmp_dict, outfile, separators=(',', ':'))

    dprint(f'puck {apid} done')


pids = list(range(1,208,2))
if __name__ == '__main__':
    start = time.time()
    with Pool(20) as p:
        p.map(process_pid, pids)
    end = time.time()
    dprint(f'Total time {end - start} seconds.')
