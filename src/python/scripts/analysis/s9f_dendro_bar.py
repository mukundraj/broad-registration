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

Supplementary:

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9f/gene_jsons_s9f gs://ml_portal2/test_data2/s9f/

References: 

https://stackoverflow.com/questions/1312331/using-a-global-dictionary-with-threads-in-python
https://stackoverflow.com/questions/38795826/optimizing-multiprocessing-pool-with-expensive-initialization

Created by Mukund on 2022-08-03

"""

import json
import sys
from produtils import dprint
import anndata as ann
import gc
from scipy.sparse import csr_matrix, hstack
import numpy as np
from allensdk.core.reference_space_cache import ReferenceSpaceCache
from multiprocessing import Pool
from multiprocessing import Manager
import functools


def initial_populate_data():
    data_root = sys.argv[1]
    ip_aggr_data_folder = data_root+sys.argv[2]
    start_pid = sys.argv[3]
    end_pid = sys.argv[4]
    op_folder = data_root+sys.argv[5]


    reference_space_key = 'annotation/ccf_2017/'
    resolution = 25
    # resolution = 10
    rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
    # rspc = ReferenceSpaceCache(resolution, reference_space_key)
    # ID 1 is the adult mouse structure graph
    allentree = rspc.get_structure_tree(structure_graph_id=1)
    name_map = allentree.get_name_map()
    # global all_region_ids
    all_region_ids = list(name_map.keys())


    # get puck ids
    pids = list(range(int(start_pid), int(end_pid)+1, 2))
    if (5 in pids):
        pids.remove(5)
    if (77 in pids):
        pids.remove(77)
    if (167 in pids):
        pids.remove(167)

    # data = {}
    processed_genes = set()
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
            processed_genes.add(gene)

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

    # exit(0)
    ## hydrate parent regions with no direct assignment of beads
    gene_idxs = list(range(len(processed_genes)))
    gene_items = zip(gene_idxs, processed_genes)
    gene_items = list(gene_items)
    return all_region_ids, gene_items, data, len(pids)

# for gene_idx, gene in enumerate(list(processed_genes)):



def process_gene(item, all_region_ids, data, len_pids):
    reference_space_key = 'annotation/ccf_2017/'
    resolution = 25
    # resolution = 10
    rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
    # rspc = ReferenceSpaceCache(resolution, reference_space_key)
    # ID 1 is the adult mouse structure graph
    allentree = rspc.get_structure_tree(structure_graph_id=1)
    gene_idx, gene= item
    data[gene].pop(str(0), None) # removing rid for beads outside tissue
    ancestors_data = {}
    # for each region id
    for ind, rid in enumerate(all_region_ids):
        if (ind==0):
            dprint("hydrating parents - gene_idx", gene_idx, ",region: ", ind, "outof", len(all_region_ids))
        ancestors_data[rid] = {"puck_dist":[int(0)] * len_pids}
        # if not leaf then loop through all possible descendents
        if not data[gene].get(rid):

            for sub_rid in data[gene].keys():
                # dprint(sub_rid,rid)
                # check if indeed descendent
                if allentree.structure_descends_from(int(sub_rid), int(rid)):
                    ances_data = np.array(ancestors_data[rid]["puck_dist"])
                    sub_data = np.array(data[gene][sub_rid]["puck_dist"])
                    # dprint(ances_data)
                    # dprint(sub_data)
                    ancestors_data[rid]["puck_dist"] = [int(x) for x in list(ances_data+sub_data)]

    # merge ancestor data dict with descendent data dict
    data[gene] = {**data[gene], **ancestors_data}

    data_root = sys.argv[1]
    op_folder = data_root+sys.argv[5]
    op_file = f'{op_folder}/{gene}.json'
    with open(op_file, 'w') as outfile:
        json.dump(data[gene], outfile, separators=(',', ':'))


data = {}
if __name__=='__main__':
    all_region_ids, gene_items, data, len_pids = initial_populate_data()
    with Pool(2) as p:
        p.map(functools.partial(process_gene, all_region_ids=all_region_ids, data=data, len_pids=len_pids), gene_items)


