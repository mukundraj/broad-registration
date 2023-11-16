"""Format aggregated single cell tab data aggregated by cellclass and clade
(generated by s2b_gen_aggr_sstab_h5ad.R) into a zarr file for use in the
portal.

Usage:

python s3a_gen_aggr_sstab_zarr.py
    in: data_root
    in: input directory
    op: output directory

Example:

python src/python/scripts/analysis_sc/s3a_gen_aggr_sstab_zarr.py \
    ~/Desktop/work/data/mouse_atlas \
    /single_cell/s2/agged_h5ad \
    /single_cell/s3 \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s3/aggedSCdata.zarr gs://bcdportaldata/batch_231112/single_cell/v3/aggedSCdata.zarr

Created by Mukund on 2023-11-15
"""
import zarr
import anndata
import numcodecs

import sys

data_root = sys.argv[1]
ipath = data_root+sys.argv[2]
op_zarr_path = data_root+sys.argv[3]+'/aggedSCdata.zarr'

fnames = [
    'clades_agged_avgs',
    'clades_agged_counts',
    'clades_agged_nz_pct',
    'cellclasses_agged_avgs',
    'cellclasses_agged_counts',
    'cellclasses_agged_nz_pct',
]

z = zarr.open( op_zarr_path, mode='w')

for fname in fnames:
    print(f'Loading {fname}')
    adata = anndata.read_h5ad(ipath+'/'+fname+'.h5ad')

    print(f'Writing {fname}')
    z[fname] = adata.X



metadataGroup = z.create_group(f'metadata', overwrite=True)

# get clade names
clade_names = []
adata = anndata.read_h5ad(ipath+'/clades_agged_avgs.h5ad')
for clade in adata.obs.to_numpy():
    clade_names.append(clade[0])

nClades = len(clade_names)
cladesArray = metadataGroup.zeros('clades', shape=(nClades), dtype='object', object_codec=numcodecs.VLenUTF8())
cladesArray[:] = clade_names

# get cell class names
cellclass_names = []
adata = anndata.read_h5ad(ipath+'/cellclasses_agged_avgs.h5ad')
for cellclass in adata.obs.to_numpy():
    cellclass_names.append(cellclass[0])
nCellclasses = len(cellclass_names)
cellclassesArray = metadataGroup.zeros('cellclasses', shape=(nCellclasses), dtype='object', object_codec=numcodecs.VLenUTF8())
cellclassesArray[:] = cellclass_names

print('Done!')
