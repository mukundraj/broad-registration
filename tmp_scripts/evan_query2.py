"""Create a csv of all celltypes and available metadata for each celltype.

Updates: 
    - 2023-12-17: new cell list based based on new scZarr file. Also removing
      additionalMetadata column.

phthon evan_query2.py
    ip: data_root
    ip: scZarr path relative to data_root
    op: output file

Example:

python tmp_scripts/evan_query2.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s1/scZarr_231207.zarr \
    /misc/for_evan/celltypes_list.csv \

Created by Mukund on 2023-10-18
"""


import sys
import zarr
from produtils import dprint
import os


data_root = sys.argv[1]
zarr_relative = sys.argv[2]
out_file = data_root+sys.argv[3]

zarr_path = data_root + zarr_relative

dprint(zarr_path)

z = zarr.open(zarr_path, mode='r')



dprint(z.tree())

cellclusters = z.obs['clusters']
cellclusters = [x.split('=')[1] for x in cellclusters]

# additionalMetadata = z.metadata['additionalMetadata']


# dprint(cellclusters[:])
# dprint(additionalMetadata[:])

# write out csv with columns cellclusters and additionalMetadata

with open(out_file, 'w') as f:
    f.write('celltype,status\n')
    for i in range(len(cellclusters)):
        f.write(cellclusters[i] + ',''\n')




