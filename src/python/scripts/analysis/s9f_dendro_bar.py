"""
Script to generate genewise jsons for dendrogram linked barplot. For each gene, 
creates a file with json containing data sorted regionwise.

Usage:

python s9f_dendro_bar.py \
    inp: data root
    inp: path to regionwise aggregated gene exp data
    inp inp: start_pid end_pid
    out: path to output folder

Usage example:

python src/python/scripts/analysis/s9f_dendro_bar.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/interim \
    1 207 \
    /data_v3_nissl_post_qc/s9_analysis/s9f/gene_jsons_s9f \

Created by Mukund on 2022-08-03

"""

import json
import sys
from produtils import dprint
import anndata as ann
import gc
from scipy.sparse import csr_matrix, hstack
import numpy as np

data_root = sys.argv[1]
ip_aggr_data_folder = data_root+sys.argv[2]
start_pid = sys.argv[3]
end_pid = sys.argv[4]
op_folder = data_root+sys.argv[5]

data = {}

# get puck ids
pids = list(range(int(start_pid), int(end_pid)+1, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

for pids_idx, pid in enumerate(pids):
    assert(pid!=5 and pid!=77 and pid!=167)
    nis_id_str = str(pid).zfill(3)

    # read puck data and get region ids
    aggr_counts_file = f'{ip_aggr_data_folder}/aggr_counts_{nis_id_str}.h5ad'
    aggr_counts = ann.read_h5ad(aggr_counts_file)


    # iterate over each gene
    # if len(genes_list)>0 and genes_list[0] !="":
    #     genes = genes_list
    # else:
    genes = list(aggr_counts.var_names)
    regions = list(aggr_counts.obs_names)

    for gene_idx, gene in enumerate(genes):

        if gene=='Pcp4' or gene=='Tph1':
            dprint(f'Found {gene} at {gene_idx}, pid: {pid}')
        else:
            continue
        if (gene_idx%500==0):
            collected = gc.collect()
            dprint('gene_idx', gene_idx, 'pid', pid, 'collected', collected)

        if gene not in data.keys():
            data[gene] = {}

        spec_gene_regagg_cnts = csr_matrix(aggr_counts.X).getcol(gene_idx)
        spec_gene_regagg_cnts_dense = np.squeeze(np.array(spec_gene_regagg_cnts.todense())).astype(int)
        if ('0' in regions):
            out_idx = regions.index('0')
            spec_gene_regagg_cnts_dense[out_idx] = 0

        obs_names_list = list(aggr_counts.obs_names)
        region_to_idx = {}
        for idx, obs in enumerate(obs_names_list):
            region_to_idx[int(obs)] = idx

        for rid in obs_names_list:
            if rid not in data[gene].keys():
                data[gene][rid] = {"puck_dist":[int(0)] * len(pids)}

            cur_gene_region_cnt = spec_gene_regagg_cnts_dense[region_to_idx[int(rid)]]
            data[gene][rid]["puck_dist"][pids_idx]=int(cur_gene_region_cnt)


# write out genewise json files after updating max_count_idx field
for gene in data.keys():
    # for rid in data[gene].keys():
    #     puck_counts = np.array(data[gene][rid]['puck_dist'])
    #     max_idx = np.argmax(puck_counts)
    #     data[gene][rid]["max_count_idx"] = int(max_idx)
    data[gene].pop(str(0), None) # removing rid for beads outside tissue
    op_file = f'{op_folder}/{gene}.json'
    with open(op_file, 'w') as outfile:
        json.dump(data[gene], outfile, separators=(',', ':'))


dprint(data)
