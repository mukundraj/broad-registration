"""
Script to generate regionwise and puckwise sorted experssion counts data for viewer.

Usage:

python s9e_sorted_counts.py \
    inp: data root
    inp: path to regionwise aggregated gene exp data
    inp: path to ccf region to region id mapping json
    inp inp: start_pid end_pid
    inp: genes list, all genes if ""
    out: path to output folder

Usage example:

python src/python/scripts/analysis/s9e_sorted_counts.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/interim \
    /data_v3_nissl_post_qc/s9_analysis/ccf_regions.json \
    1 3 \
    "Pcp4, Gad2" \
    /data_v3_nissl_post_qc/s9_analysis/s9e/gene_jsons_s9e \

Created by Mukund on 2022-07-12

Supplementary:

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9e/gene_jsons_s9e gs://ml_portal2/test_data2/

"""

import sys
from produtils import dprint
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix, hstack
import json
import gc

data_root = sys.argv[1]
ip_aggr_data_folder = data_root+sys.argv[2]
ccf_name_to_id_json = data_root+sys.argv[3]
start_pid = sys.argv[4]
end_pid = sys.argv[5]
genes_list_str = sys.argv[6]
op_folder = data_root+sys.argv[7]

# get gene list
genes_list = genes_list_str.split(",")

# get puck ids
pids = list(range(int(start_pid), int(end_pid)+1, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

# get region ids


region_aggred_counts = {}
puck_aggred_counts = {}
gene_region_ids = {}

# iterate over each puck
for pid in pids:
    assert(pid!=5 and pid!=77 and pid!=167)
    nis_id_str = str(pid).zfill(3)

    # read puck data and get region ids
    aggr_counts_file = f'{ip_aggr_data_folder}/aggr_counts_{nis_id_str}.h5ad'
    aggr_counts = ann.read_h5ad(aggr_counts_file)


    # iterate over each gene
    if len(genes_list)>0 and genes_list[0] !="":
        genes = genes_list
    else:
        genes = list(aggr_counts.var_names)

    for gene_idx, gene in enumerate(genes):
        # if gene=='Pcp4' or gene=='Gad2':
        # if (len(genes_list)>0) and gene in genes_list:
        #     dprint(f'Found {gene} at {gene_idx}')
        # else:
        #     continue
        if (gene_idx%500==0):
            collected = gc.collect()
            dprint('gene_idx', gene_idx, 'pid', pid, 'collected', collected)

        if gene not in region_aggred_counts.keys():
            region_aggred_counts[gene] = {}
            gene_region_ids[gene] = set()
            puck_aggred_counts[gene] = {}

        spec_gene_regagg_cnts = csr_matrix(aggr_counts.X).getcol(gene_idx)
        spec_gene_regagg_cnts_dense = np.squeeze(np.array(spec_gene_regagg_cnts.todense())).astype(int)
        # dprint(np.shape(spec_gene_regagg_cnts_dense))

        obs_names_list = list(aggr_counts.obs_names)
        region_to_idx = {}
        for idx, obs in enumerate(obs_names_list):
            region_to_idx[int(obs)] = idx

        # iterate over each region in puck
        for rid in obs_names_list:
            gene_region_ids[gene].add(int(rid))
            rkey = (pid, int(rid))

            cur_gene_region_cnt = spec_gene_regagg_cnts_dense[region_to_idx[int(rid)]]
            # if key not exist, add key, else increment count
            if (rkey not in region_aggred_counts[gene].keys()):
                region_aggred_counts[gene][rkey] = cur_gene_region_cnt
            else:
                region_aggred_counts[gene][rkey] += cur_gene_region_cnt

            pkey = (pid, -1) # -1 region added to maintain structure as region list
            if (pkey not in puck_aggred_counts[gene].keys()):
                puck_aggred_counts[gene][pkey] = cur_gene_region_cnt
            else:
                puck_aggred_counts[gene][pkey] += cur_gene_region_cnt


# create a map from region ids to region names
with open(ccf_name_to_id_json ) as json_file:
    ccf_name_to_id = json.load(json_file)

region_id_to_name = {v: k for k, v in ccf_name_to_id.items()}

# write out aggregated count info per gene
for gene in region_aggred_counts.keys():

    # get and sort values in region_aggred_counts
    reg_aggred_vals = [{"key":key, "cnt": int(region_aggred_counts[gene][key])} for key in region_aggred_counts[gene].keys() if region_aggred_counts[gene][key]>0]
    reg_aggred_vals = sorted(reg_aggred_vals, key=lambda x: x['cnt'], reverse=True)

    puck_aggred_vals = [{"key":key, "cnt": int(puck_aggred_counts[gene][key])} for key in puck_aggred_counts[gene].keys()]
    puck_aggred_vals = sorted(puck_aggred_vals, key=lambda x: x['key'][0])

    gene_region_ids_to_name = [{'rid':k,'name':region_id_to_name[k]} for k in region_id_to_name.keys() if k in gene_region_ids[gene]]
    # dprint(gene_region_ids_to_name)
    # get and sort values in puck_aggred_counts
    out_data = {"sorted_regionwise_cnts":reg_aggred_vals, "sorted_puckwise_cnts":puck_aggred_vals, "regionid_to_name":gene_region_ids_to_name}

    # write out outdata as json
    op_file = f'{op_folder}/{gene}.json'
    with open(op_file, 'w') as outfile:
        json.dump(out_data, outfile, separators=(',', ':'))

dprint("done")
