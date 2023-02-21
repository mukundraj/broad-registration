""" Creates folder with anndata files with avg and percent values and a readme
file as well. Also, will contain all metadata in same anndata files incl bead
coords, gene expression across beads as gene metadata.

Caveat - would neet about 50GB space for output. Alternatively, can move files
to Gdrive as being created. Each file creation can take over a minute providing
plenty time to move to Gdrive.

Usage:

python s1c_singlecell.py
    io: data root
    ip: scZarr file path with avg and percent values
    op: op_folder

Example:

python src/python/scripts/misc/for_downloads/s1c_singlecell.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /single_cell/s1/scZarr.zarr \
    /misc/for_downloads/s1c_single_cell \

Supplemantary:

copy op_folder/scdata/* to ...Gdrive/braincelldata.org/data_single_cell_metadata

Created by Mukund on 2023-01-25
"""

from produtils import dprint
import anndata
import sys
import os
import shutil
import zarr
import numpy as np
from scipy.sparse import csr_matrix



from pathlib import Path
path_root = Path(__file__).parents[5]
sys.path.append(str(path_root))
import src.python.utils.filename_maps as fmaps

data_root = sys.argv[1]
scZarr_file = f'{data_root}/{sys.argv[2]}'
op_folder = f'{data_root}/{sys.argv[3]}'

pidToSrno, srnoToPid = fmaps.get_pidToSrno()

# create a dicrectory scdata in opfolder and overwrite if exists
scdata_folder = f'{op_folder}/singlecell_data'
dprint(scdata_folder)
if os.path.exists(scdata_folder):
    shutil.rmtree(scdata_folder)
os.mkdir(scdata_folder)

dprint(scZarr_file)
# read scZarr file
scZarr = zarr.open(scZarr_file, mode='r')

print(scZarr.tree())

# copy readmen.txt from s1d_single_cell folder to scdata folder
readme_file = f'{op_folder}/readme.txt'
shutil.copy(readme_file, scdata_folder)

for srno in range(1, 102):
# for srno in range(1, 3):
    pid = srnoToPid[str(srno).zfill(3)]

    adfile = f'{scdata_folder}/sr_{srno}.h5ad'

    # dprint(srno, pid, adfile)

    # read avg and percent values from zarr
    avg = csr_matrix(scZarr['avg'].X[:], dtype=np.float32)
    pct = csr_matrix(scZarr['nz_pct'].X[:], dtype=np.float32)
    dprint(type(avg))
    dprint(avg.shape)

    # add cell metadata to anndata file
    an_fname = f'{scdata_folder}/sr_{srno}.h5ad'

    # initialize the anndata file with avg values and pct values
    an = anndata.AnnData(X=avg)
    an.layers['pct'] = pct


    an.obs['topstructs'] = scZarr['metadata'].topstructs
    an.obs['cellclass'] = scZarr['metadata'].cellclasses
    # add cell scores as cell metadata to anndata file
    # an.obsm['cell_scores'] = np.zeros((5030, 3))
    # add bead index for mapped celltypes

    # add gene counts across beads as gene metadata to anndata file

    # add the dissectate info metadata to anndata file (as a dict)
    # # an.uns["dissectates"] = {'cellname': ['D1', 'D2']}
    # an.uns["bead_coords"] =  np.zeros((20000,3), dtype=np.int16)

    # add cellclass metadata to anndata file

    # display overview
    dprint(an)

    # write out an_fname
    an.write(an_fname, compression="gzip")



