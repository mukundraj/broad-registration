""" Creates a mapping from slide number to names of CCF (leaf) regions in that
slice.

Usage:

python src/python/scripts/misc/for_alyssa/gen_slide2ccf.py \
    ~/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9g/nz_aggr_counts_cshl_230222 \
    > ~/Desktop/slide_to_ccfregion_mapping.tsv

Generate a table that maps slide to CCF region name.

Created by Mukund on 2024-05-14
"""

from allensdk.core.reference_space_cache import ReferenceSpaceCache
import sys
import os

from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import src.python.utils.filename_maps as fmaps

# path_to_csv_with_regioninfo = '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9g/nz_aggr_counts_cshl_230222'

path_to_csv_with_regioninfo = sys.argv[1]

# for name map

reference_space_key = 'annotation/ccf_2017'
resolution = 25
rspc = ReferenceSpaceCache(resolution, reference_space_key, manifest='manifest.json')
# ID 1 is the adult mouse structure graph
tree = rspc.get_structure_tree(structure_graph_id=1) 

id_to_name_map = tree.get_name_map()

# print(id_to_name_map)

pidToSrno, srnoToPid = fmaps.get_pidToSrno()

# print (pidToSrno)

# get all csv file names in path_to_csv_with_regioninfo

files = os.listdir(path_to_csv_with_regioninfo)
csv_files = [f for f in files if f.endswith('.csv')]
# print(csv_files)

# sort the csv files
csv_files = sorted(csv_files)

# print header of the table
print('sno\tregion_name')

# iterate over each csv file and read the table
for csv_file in csv_files:
    # print (csv_file)
    # get pid
    pid = int(csv_file.split('_')[-1].split('.')[0])
    # print (pid)
    with open(f'{path_to_csv_with_regioninfo}/{csv_file}', 'r') as f:
        lines = f.readlines()
        # skip the first line
        lines = lines[1:]
        for line in lines:
            line = line.strip().split(',')
            # print (line)
            regionid = int(line[0])
            count = int(line[1])
            if (regionid == 0):
                continue
            region_name = id_to_name_map[regionid]
            srno = pidToSrno[pid]
            print(f'{srno}\t{region_name}')

