"""
Script to cellrate regionwise and puckwise sorted scores data for
viewer's CellSpatial tab's interactive barplot.

Usage:

python s2c_sorted_counts.py \
    i/o: data_root \
    inp: path to regionwise aggregated cell exp data
    inp: path to ccf region to region id mapping json
    inp inp: start_pid end_pid
    out: path to output folder

Usage example:

python src/python/scripts/analysis_sc/s2c_sorted_counts.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s1/s1d_region_agg \
    /data_v3_nissl_post_qc/s9_analysis/ccf_regions.json \
    1 207 \
    /cell_spatial/s2/s2c/cell_jsons_s2c \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/cell_spatial/s2/s2c/cell_jsons_s2c gs://bcdportaldata/cellspatial_data/freqbars/cell_jsons_s2c

Created by Mukund on 2022-10-27

"""


import sys
from produtils import dprint
import anndata as ann
import numpy as np
from scipy.sparse import csr_matrix, hstack
import json
import gc
from pathlib import Path
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.allen as allen

data_root = sys.argv[1]
ip_aggr_data_folder = data_root+sys.argv[2]
ccf_name_to_id_json = data_root+sys.argv[3]
start_pid = sys.argv[4]
end_pid = sys.argv[5]
# genes_list_str = sys.argv[6]
op_folder = data_root+sys.argv[6]

# get gene list
# genes_list = genes_list_str.split(",")

# get puck ids
pids = list(range(int(start_pid), int(end_pid)+1, 2))
if (5 in pids):
    pids.remove(5)
if (77 in pids):
    pids.remove(77)
if (167 in pids):
    pids.remove(167)

# get region ids
agg_rids, agg_rnames = allen.get_allen_ish_viewer_regions()


region_aggred_counts = {}
puck_aggred_counts = {}
cell_region_ids = {}



# iterate over each puck
for pid in pids:
    assert(pid!=5 and pid!=77 and pid!=167)
    nis_id_str = str(pid).zfill(3)

    # read puck data and get region ids
    aggr_counts_file = f'{ip_aggr_data_folder}/aggr_scores_{nis_id_str}.h5ad'
    aggr_counts = ann.read_h5ad(aggr_counts_file)


    # iterate over each gene
    # if len(genes_list)>0 and genes_list[0] !="":
    #     genes = genes_list
    # else:
    cells = list(aggr_counts.var_names)
    regions = list(aggr_counts.obs_names)
    csr_cell_reagg_cnts = csr_matrix(aggr_counts.X)
    # all_cell_count = csr_cell_reagg_cnts.sum()
    # normalizer_val = csr_cell_reagg_cnts.sum()/100000 # to get counts per 10K
    normalizer_val = 1

    for cell_idx, cell in enumerate(cells):

        # if cell=='Pcp4' or cell=='Tph1':
        #     dprint(f'Found {cell} at {cell_idx}, pid: {pid}')
        # else:
        #     continue
        if (cell_idx%500==0):
            collected = gc.collect()
            dprint('cell_idx', cell_idx, 'pid', pid, 'collected', collected)

        if cell not in region_aggred_counts.keys():
            region_aggred_counts[cell] = {}
            cell_region_ids[cell] = set()
            puck_aggred_counts[cell] = {}


        spec_cell_regagg_cnts = csr_cell_reagg_cnts.getcol(cell_idx)
        # spec_cell_regagg_cnts_dense = np.squeeze(np.array(spec_cell_regagg_cnts.todense())).astype(int)
        spec_cell_regagg_cnts_dense = np.squeeze(np.array(spec_cell_regagg_cnts.todense())).astype(float)
        # if (cell.split("=")[1]=='Inh_Pax6_Nkx2-2_1'):
        #     dprint('cell_idx', cell_idx, 'pid', pid, 'spec_cell_regagg_cnts_dense', spec_cell_regagg_cnts_dense, 'shape', spec_cell_regagg_cnts_dense.shape)
        # dprint(np.shape(spec_cell_regagg_cnts_dense))
        if ('0' in regions):
            out_idx = regions.index('0')
            # dprint('out agg: ', spec_cell_regagg_cnts_dense[out_idx])
            spec_cell_regagg_cnts_dense[out_idx] = 0


        obs_names_list = list(aggr_counts.obs_names)
        region_to_idx = {}
        for idx, obs in enumerate(obs_names_list):
            region_to_idx[int(obs)] = idx

        # dprint(obs_names_list)
        # iterate over each region in puck
        for rid in obs_names_list:
            cell_region_ids[cell].add(int(rid))

            cur_cell_region_cnt = spec_cell_regagg_cnts_dense[region_to_idx[int(rid)]]

            # if key not exist, add key, else increment count
            # rkey = (pid, int(rid))
            # if (rkey not in region_aggred_counts[cell].keys()):
            #     region_aggred_counts[cell][rkey] = cur_cell_region_cnt
            # else:
            #     region_aggred_counts[cell][rkey] += cur_cell_region_cnt

            pkey = (pid, -1) # -1 region added to maintain structure as region list
            # dprint(pkey, rid, cur_cell_region_cnt)
            if (pkey not in puck_aggred_counts[cell].keys()):
                puck_aggred_counts[cell][pkey] = cur_cell_region_cnt / normalizer_val
            else:
                puck_aggred_counts[cell][pkey] += cur_cell_region_cnt / normalizer_val

            # iterate over each agg_region - if rid in aggregion, increment aggregion cnt
            for agridx, aggrids in agg_rids.items():
                if int(rid) in aggrids:
                    rkey = (pid, agridx)
                    if (rkey not in region_aggred_counts[cell].keys()):
                        region_aggred_counts[cell][rkey] = cur_cell_region_cnt / normalizer_val
                    else:
                        region_aggred_counts[cell][rkey] += cur_cell_region_cnt / normalizer_val



# reduce region_aggred_counts across pucks into only one puck (with most count)

for cell in region_aggred_counts.keys():
    region_vals = list([{'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0},
                        {'max':-1, 'maxpid':-1, 'sum':0}])

    for key,val in region_aggred_counts[cell].items():
        region_vals[key[1]]['sum'] += val
        if val > region_vals[key[1]]['max']:
            region_vals[key[1]]['max'] = val
            region_vals[key[1]]['maxpid'] = key[0]

    region_aggred_counts[cell] = {}
    # region_vals[0] = {}
    for agridx,item in enumerate(region_vals):
        rkey = (item['maxpid'],agridx)
        region_aggred_counts[cell][rkey] = item['max']


# dprint(region_aggred_counts)

# exit(0)

# create a map from region ids to region names
with open(ccf_name_to_id_json ) as json_file:
    ccf_name_to_id = json.load(json_file)

region_id_to_name = {v: k for k, v in ccf_name_to_id.items()}

# write out aggregated count info per cell
for cell in region_aggred_counts.keys():


    puck_aggred_vals = [{"key":key, "cnt": round(float(puck_aggred_counts[cell][key]), 4)} for key in puck_aggred_counts[cell].keys()]
    puck_aggred_vals = sorted(puck_aggred_vals, key=lambda x: x['key'][0])
    puck_aggred_vals = [{**item, 'sr':sr} for sr,item in enumerate(puck_aggred_vals)]

    pid_to_sr = {}
    for item in puck_aggred_vals:
        pid_to_sr[item['key'][0]] = item['sr']

    # get and sort values in region_aggred_counts
    reg_aggred_vals = [{"key":key, "cnt": round(float(region_aggred_counts[cell][key]), 4)} for key in region_aggred_counts[cell].keys() if region_aggred_counts[cell][key]>0]
    # reg_aggred_vals = sorted(reg_aggred_vals, key=lambda x: x['cnt'], reverse=True)
    reg_aggred_vals = [{**item, 'nm':agg_rnames[i], 'sr': pid_to_sr[item['key'][0]]} for i, item in enumerate(reg_aggred_vals)]

    # cell_region_ids_to_name = [{'rid':k,'name':region_id_to_name[k]} for k in region_id_to_name.keys() if k in cell_region_ids[cell]]

    out_data = {"regionwise_cnts":reg_aggred_vals, "sorted_puckwise_cnts":puck_aggred_vals}

    # write out outdata as json
    cellname = cell.split('=')[1]
    op_file = f'{op_folder}/{cellname}.json'
    with open(op_file, 'w') as outfile:
        
        json.dump(out_data, outfile, separators=(',', ':'))

dprint("done")
