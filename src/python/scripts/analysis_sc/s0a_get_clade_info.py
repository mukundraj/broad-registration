"""
Generates a file with formatted clade info.

Usage:

python s0a_get_clade_info.py \
    ip: data root
    ip: raw clade file
    ip: file containing clusterNames and sizes
    op: processed clade file
    op: list of cells not mapped to any clade

Example:

python src/python/scripts/analysis_sc/s0a_get_clade_info.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s0/CellType_Clades.csv \
    /single_cell/s0/raw_v2/20220912_QC_summary/clusterSize.csv \
    /single_cell/s1/processed_clade_info.csv \
    /single_cell/s1/unmapped_cells.csv \

Created by Mukund on 2023-10-04
"""

import sys
import csv
from produtils import dprint

data_root = sys.argv[1]
ip_clade_file = f'{data_root}/{sys.argv[2]}'
clustersize_csv_file = data_root+sys.argv[3]
op_clade_file = f'{data_root}{sys.argv[4]}'
op_unmapped_cells = f'{data_root}{sys.argv[5]}'

clusterNames = []

# first, read clusterSize file
with open(clustersize_csv_file, 'r') as f:
    # pass the file object to reader() to get the reader object
    csv_reader = csv.reader(f)
    # Iterate over each row in the csv using reader object
    for idx,row in enumerate(csv_reader):
        # row variable is a list that represents a row in csv
        if (idx>0):
            clusterName = row[0].split('=')[1]
            clusterNames.append(clusterName)

dprint(len(clusterNames))

# read clade file
cell_to_clade_map = {}
with open(ip_clade_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    next(reader) # skip header
    for row in reader:
        # dprint(row[0], row[2], row[7], row[8])
        leafs = row[2].split(sep='|')
        # dprint(row)
        for leaf in leafs:
            cell_to_clade_map[leaf] = [row[0],row[1]]

unmapped_cells = set()
dprint(len(cell_to_clade_map.keys()))

for cell in clusterNames:
    if (cell in cell_to_clade_map.keys()):
        continue
    else:
        unmapped_cells.add(cell)


# write out processed clade file
with open(op_clade_file, 'w') as f:
    writer = csv.writer(f)
    for cell in clusterNames:
        if (cell in cell_to_clade_map.keys()):
            writer.writerow([cell, cell_to_clade_map[cell][0], cell_to_clade_map[cell][1]])
        else:
            writer.writerow([cell, 'Clade0', ''])


dprint(len(unmapped_cells))

# write out list of cells that are missing a clade
with open(op_unmapped_cells, 'w') as f:
    writer = csv.writer(f)
    for umcell in unmapped_cells:
        writer.writerow([umcell])

