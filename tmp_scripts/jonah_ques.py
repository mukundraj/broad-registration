"""
Misc code to answer Jonah's query. Get total number of beads with score > 0.3
in for gene  Inh_Frmd7_Lamp5

python tmp_scripts/jonah_ques.py \

Created by Mukund on 2023-09-21 
"""


import zarr
import json
import numpy as np


pids = list(range(1,208,2))

total = 0
for pid in pids:

    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    print(f'starting apid {apid}..................')

    puckfolder = f'/Users/mraj/Desktop/work/data/mouse_atlas/cell_spatial/s1/cellspatial_data/cellscores_cshl/puck{apid}'

    zarrfile = f'{puckfolder}/cellxbead.zarr'

    optionsfile = f'{puckfolder}/cellOptions.json'

    # print(zarrfile)
    # print(optionsfile)

    z = zarr.open(zarrfile)
    # print(z.tree())

    f = open(optionsfile)
    cells = json.load(f)['cellOptions']
    f.close()

    # print(len(cells))

    celltype = 'Inh_Frmd7_Lamp5'

    if celltype in cells:
        index = cells.index(celltype)
        # print(index)

        celldata = np.array(z.X[index,:])
        # print(np.sum(celldata))
        celldata[celldata<0.3] = 0

        nonzero = np.count_nonzero((celldata))
        # print(nonzero)
        total += nonzero


    print ('tillpid', apid, 'total', total)





