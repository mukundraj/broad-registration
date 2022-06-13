"""

Usage:

Usage example:

python src/python/scripts/analysis/s9c_qc_main.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s9_analysis/aggregated_labels \
    /Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats \
    /temp \
    /data_v3_nissl_post_qc/s9_analysis/qc_mats/processed \
    "Ndnf" \
    79 83 \
    False \
    /data_v3_nissl_post_qc/s9_analysis/qc_mats/processed

Created by Mukund on 2022-06-11

"""

import sys
import subprocess
from produtils import dprint
from allensdk.core.reference_space_cache import ReferenceSpaceCache
import os
import json
import numpy as np
from pathlib import Path
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.allen as allen

data_root = sys.argv[1]
ip_folder_labels = data_root+sys.argv[2]
ip_folder_counts = sys.argv[3]
op_folder_tmp = data_root+sys.argv[4]
op_folder = data_root+sys.argv[5]
genes_list_str = sys.argv[6]
start_pid = sys.argv[7]
end_pid = sys.argv[8]
nonzero = sys.argv[9]=='True' # note its == and not =
op_folder = sys.argv[1]+sys.argv[10]

def stringToList(string):
    listRes = list(string.split(" "))
    return listRes

genes_list = stringToList(genes_list_str)



# regions_names = [
#    [ "Primary somatosensory area, trunk, layer 1"],
#    [ "Primary somatosensory area, trunk, layer 2/3"],
#    [ "Primary somatosensory area, trunk, layer 4"],
#    [ "Primary somatosensory area, trunk, layer 5"],
#    [ "Primary somatosensory area, trunk, layer 6a"],
#    [ "Primary somatosensory area, trunk, layer 6b"],
#    [ "Hippocampal formation"]]

# reference_space_key = 'annotation/ccf_2017'
# resolution = 25
# rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# # ID 1 is the adult mouse structure graph
# tree = rspc.get_structure_tree(structure_graph_id=1)

# regions_ids = []
# for layer_names in regions_names:
#     dprint(layer_names)
#     # dprint(tree.get_structures_by_name(layer_names))
#     regions_ids.append(list(map(lambda x: x['id'], tree.get_structures_by_name(layer_names))))


# regions_ids = allen.get_cortex_layer_ids_lists()

# regions_ids_str = ''
# for layer in regions_ids:
#     for id in layer:
#         regions_ids_str+=str(id)+','
#     regions_ids_str = regions_ids_str[:-1]
#     regions_ids_str+=';'
# regions_ids_str = regions_ids_str[:-1]
# dprint(regions_ids_str)
# regions_ids_str = '68,107;667,219;12995,865;229,648;644,947'

regions_ids_str = '1;2;3;4;5'

# genes_list = ['Dcn', 'Ndnf']
# genes_list = ['Gad1']
# genes_list = ['Ndnf']
# region_ids = '234,234,234;33,33,32;'

# subprocess.run(["fd", "^.*tif$" , "-x", "convert", "{}", "-alpha", "off", "{}" ], cwd=op_path)


# iterate over genes_list

for gene in genes_list:

    # call R script
    r_op_file_name = f'{op_folder_tmp}/{gene}.csv'
    dprint(ip_folder_labels)
    dprint(ip_folder_counts)
    dprint(r_op_file_name)
    subprocess.run(["Rscript", "src/python/scripts/analysis/s9c_qc_prep.R", \
                    gene, ip_folder_labels, ip_folder_counts, r_op_file_name, \
                    regions_ids_str, \
                    start_pid, end_pid, \
                    str(nonzero)])

    # read gene matrix and prepare dict object
    data_z = np.genfromtxt(r_op_file_name, delimiter=',')
    n_regions,n_pucks = np.shape(data_z)
    # dprint(data_z.tolist())
    dprint(n_regions,n_pucks)
    data_dict = {"x": list(range(1, n_pucks+1)), "y": list(range(1, n_regions+1)), "z":data_z.tolist()}

    # write out gene data json
    fname_val = "fnz" if nonzero else "fumi"
    fname = f'{fname_val}_{gene}.json'
    op_gene_data_file = f'{op_folder}/{fname}'
    dprint(nonzero, fname_val)
    with open(op_gene_data_file, 'w') as outfile:
        json.dump(data_dict, outfile)

    # prepare metadata
    desc_val = "nonzero count" if nonzero else "count"
    metadata_dict = {"name":f'{gene}',
                     "desc":f'gene {desc_val} normalized by nUMI in region',
                     "filename":f'{fname}'}

    dprint(metadata_dict)

    # update metadata file
    metadata_path =  f'{op_folder}/metadata.json'
    file_exists = os.path.exists(metadata_path)

    meta_key_val = 'ana_metadata_nonzero' if nonzero else 'analysis_metadata'
    if (file_exists):
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        dprint(metadata[meta_key_val])
        idxs = list(i for i, x in enumerate(metadata[meta_key_val]) if x['name'] == 'Ndnf')
        if (len(idxs)>0):
            metadata[meta_key_val][idxs[0]] = metadata_dict
        else:
            metadata[meta_key_val].append(metadata_dict)
    else:
        metadata = {
            "analysis_metadata": [],
            "ana_metadata_nonzero": []
        }
        metadata[meta_key_val].append(metadata_dict)

    # write out metadata
    with open(metadata_path, 'w') as outfile:
        json.dump(metadata, outfile)

    # copy folder to gcp storage
    # gsutil -m cp -r /Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/qc_mats/processed gs://ml_portal2/test_data2/qc_mats/

# https://github.com/conda/conda/issues/9589
# https://discuss.pytorch.org/t/intel-mkl-fatal-error-cannot-load-libmkl-core-dylib/67710/3s


