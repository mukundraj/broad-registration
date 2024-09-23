"""
Determine region for each coord and append to coords file. Creates output file
with _region.csv suffix in same folder as input coords file.

python  03_get_region.py \
    in: coords file with x, y coordinates in chuck space
    in: folder with boundary geojsons (in chuck space)
    in: nissl slice number from braincelldata.org to get regions info from 
     

python  src/python/scripts/misc/for_jonahs/03_get_region.py \
     "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2/tfmed.csv" \
     "/Users/mraj/Desktop/work/data/mouse_atlas/v3/s2/geojsons_tfmed_chksp" \
     48 \
    

Created by Mukund on 2024-09-20
"""

import sys
from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import os
import subprocess

import src.python.utils.filename_maps as fm
import geopandas
import numpy as np


tfmed_coords = sys.argv[1]

jsons_folder = sys.argv[2]

# slice number of nissl
sno = int(sys.argv[3])

# read coords
coords = np.loadtxt(tfmed_coords, delimiter=',')



# id coords of interest

pidToSrno, srnoToPid = fm.get_pidToSrno()
sno = str(sno).zfill(3)
pid = srnoToPid[sno]
print(pid)

bounds_file = f'{jsons_folder}/pid_{pid}.geojson'
bounds = geopandas.read_file(bounds_file)


print(bounds)

regions = []

total = len(coords)

print('coords')
# coords = np.array([[2080, 3606-2234],[2489, 3606-1500], [2520, 1520]])
print(coords)

# subtract y val from 3606 for all coords
coords = np.array(coords)
coords = np.array([coords[:, 0], 3606 - coords[:, 1]]).T


# iterate over points in coords
for i, x in enumerate(coords):
    region_found = False
    # check if point x is in any of the polygons 
    for j, row in bounds.iterrows():
        if row.geometry != None:
            if row.geometry.contains(geopandas.points_from_xy([x[0]], [x[1]])[0]):
                # print(row['name'])

                regions.append(row['name'])
                region_found = True
                break
    if region_found == False:
        regions.append('NA')

    print (i, "out of", total, regions[-1])

# write out coords with new column region
print(regions)
out_file = tfmed_coords.replace('.csv', f'_with_regions_{sno}.csv')
np.savetxt(out_file, np.column_stack((coords, regions)), delimiter='|', fmt='%s')

# x = [1000,1000]








