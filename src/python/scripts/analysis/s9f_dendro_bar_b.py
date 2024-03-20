"""
Second script in series to create dendro bars. Needs to be run after s9f_dendro_bar_a.py but before s9f_dendro_bar_b.py. Previously run with 240GB ram and 240GB swap space.

Usage example:

python src/python/scripts/analysis/s9f_dendro_bar_b.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9d/interim_cshl \
    1 207 \
    /data_v3_nissl_post_qc/s9_analysis/s9f/gene_jsons_s9f \

Supplementary:

gsutil -m cp -r ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9f/gene_jsons_s9f gs://ml_portal2/test_data2/s9f/

Created by Mukund on 2022-08-21

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
import pickle

data_root = sys.argv[1]
ip_aggr_data_folder = data_root+sys.argv[2]
start_pid = sys.argv[3]
end_pid = sys.argv[4]
op_folder = data_root+sys.argv[5]

with open(op_folder+'/all_region_ids.obj', 'rb') as f:
        all_region_ids = pickle.load(f);
        dprint("read all_region_ids")
with open(op_folder+'/gene_items.obj', 'rb') as f:
    gene_items = pickle.load(f);
    dprint("read gene_items")
with open(op_folder+'/data.obj', 'rb') as f:
    data = pickle.load(f);
    dprint("read data")
with open(op_folder+'/len_pids.obj', 'rb') as f:
    len_pids = pickle.load(f);
    dprint("read len_pids")

dprint(gene_items)
processed_genes = [x for (i,x) in gene_items]
dprint(processed_genes)
reference_space_key = 'annotation/ccf_2017/'
resolution = 25
# resolution = 10
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# rspc = ReferenceSpaceCache(resolution, reference_space_key)
# ID 1 is the adult mouse structure graph
allentree = rspc.get_structure_tree(structure_graph_id=1)

pids = list(range(int(start_pid), int(end_pid)+1, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

total_genes = len(processed_genes)

for gene_idx, gene in enumerate(list(processed_genes)):
    data[gene].pop(str(0), None) # removing rid for beads outside tissue
    ancestors_data = {}
    # for each region id
    for ind, rid in enumerate(all_region_ids):
        if (ind%1000==0):
            dprint("hydrating parents - gene_idx", gene_idx, "/", total_genes,",region: ", ind, "outof", len(all_region_ids))
            sys.stdout.flush()
        ancestors_data[rid] = {"puck_dist":[int(0)] * len(pids)}
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
                    ancestors_data[rid]["puck_dist"] = [round(float(x),4) for x in list(ances_data+sub_data)]

    # merge ancestor data dict with descendent data dict by adding puck_dist
    # data[gene] = {**data[gene], **ancestors_data}
            # data[gene] = {**data[gene], **ancestors_data}
            for rid in ancestors_data.keys():
                str_rid = str(rid)
                data[gene][str_rid] = ancestors_data[rid]


# write out data dict to pickle file
with open(op_folder+'/data_hydrated.obj', 'wb') as f:
    pickle.dump(data, f);
    dprint("wrote data_hydrated.obj")


# # convert all puck_dist values to string to rounding issue in float array in viewer
# for gene in data.keys():
#     for rid in data[gene].keys():
#         data[gene][rid]["puck_dist"] = [str(x) for x in data[gene][rid]["puck_dist"]]

# # write out genewise json files after updating max_count_idx field
# for gene in data.keys():
#     # for rid in data[gene].keys():
#     #     puck_counts = np.array(data[gene][rid]['puck_dist'])
#     #     max_idx = np.argmax(puck_counts)
#     #     data[gene][rid]["max_count_idx"] = int(max_idx)
#     op_file = f'{op_folder}/{gene}.json'
#     with open(op_file, 'w') as outfile:
#         json.dump(data[gene], outfile, separators=(',', ':'))
