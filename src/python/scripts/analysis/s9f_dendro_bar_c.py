"""

Usage example:

python src/python/scripts/analysis/s9f_dendro_bar_c.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/s9f \
    /data_v3_nissl_post_qc/s9_analysis/s9f/gene_jsons_s9f \

Created by Mukund on 2024-03-14
"""
import sys
import pickle
from produtils import dprint
import json

data_root = sys.argv[1]
ip_folder = data_root+sys.argv[2]
op_folder = data_root+sys.argv[3]

# hydrated data pickle file
data_hydrated_file = ip_folder+'/data_hydrated.obj'

# read hydrated data
with open(data_hydrated_file, 'rb') as f:
    data = pickle.load(f);
    dprint("read data_hydrated")

genes = list(data.keys())
# genes = genes[23255:]
# convert all puck_dist values to string to rounding issue in float array in viewer
for idx, gene in enumerate(genes):
    dprint(f"processing gene {gene} {idx+1}/{len(data.keys())}")
    current_gene_data = data[gene]
    # for rid in data[gene].keys():
    #     data[gene][rid]["puck_dist"] = [str(round(x,3)) for x in data[gene][rid]["puck_dist"]]
    current_gene_data_rids = list(current_gene_data.keys())

    current_gene_data_tmp = {}
    for rid in current_gene_data_rids:
        current_gene_data_tmp[rid] = {}

    for rid in current_gene_data_rids:
        current_gene_data_tmp[rid]["puck_dist"] = [str(round(x,4)) for x in current_gene_data[rid]["puck_dist"]]

    current_gene_data = current_gene_data_tmp


    dprint(f"writing gene {gene} {idx+1}/{len(data.keys())}", "num rids", len(current_gene_data_rids))
    op_file = f'{op_folder}/{gene}.json'
    with open(op_file, 'w') as outfile:
        json.dump(current_gene_data, outfile, separators=(',', ':'))

# # write out genewise json files after updating max_count_idx field
# for idx, gene in enumerate(data.keys()):
#     # for rid in data[gene].keys():
#     #     puck_counts = np.array(data[gene][rid]['puck_dist'])
#     #     max_idx = np.argmax(puck_counts)
#     #     data[gene][rid]["max_count_idx"] = int(max_idx)
#     dprint(f"writing gene {gene} {idx+1}/{len(data.keys())}")
#     op_file = f'{op_folder}/{gene}.json'
#     with open(op_file, 'w') as outfile:
#         json.dump(data[gene], outfile, separators=(',', ':'))
