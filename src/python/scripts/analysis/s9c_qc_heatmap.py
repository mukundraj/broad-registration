"""
Script for generating QC heatmap data with ratios involving total counts and percent non zero counts.

Usage:

python s9c_qc_heatmap.py
    inp: path to data/annex root \
    out: relative path to output dir


Usage example:

python src/python/scripts/analysis/s9c_qc_heatmap.py \
    /Users/mraj/Desktop/work/data/mlab-annex/ \
    mouse_atlas/data_v3_nissl_post_qc/s9_analysis/qc_mats/

Supplementary:
- gsutil -m cp -r ~/Desktop/work/data/mlab-annex/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/qc_mats gs://ml_portal2/test_data2/

"""
from produtils import dprint
import sys
import json


data_root = sys.argv[1]
out_path = data_root+sys.argv[2]

meta_data = {"analysis_metadata":[{"name":"sample1", "filename":"sample1.json"},
                                  {"name":"sample2", "filename":"sample2.json"}]
            }

heat_map_data = [{"data": {"x":[1,2,3], "y":[1,2], "z":[[1, 2, 3], [4, 5, 6]]}},
                 {"data": {"x":[1,2], "y":[1,2], "z":[[1, 2], [1,2]]}}]

with open(out_path+"metadata.json", "w") as outfile:
    json.dump(meta_data, outfile)


meta_data_array = meta_data["analysis_metadata"]
for idx, meta_data_item in enumerate(meta_data_array):

    with open(out_path+meta_data_item['filename'], "w") as outfile:
        json.dump(heat_map_data[idx]['data'], outfile)

