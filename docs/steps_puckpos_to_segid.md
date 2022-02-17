
## Steps for processing slide seq data, nissl data, and obtaining region id corresponding to puck positions

- **Step** : use following script to downsample (using "-series" option in bfconvert) and convert to png/tiff - `./src/python/scripts/convert_from_vsi.py`

- **Step** : generate slide seq data using following script - `./src/scripts/process-rds-files.R`

- **Step** : downsample slide seq data using following script - `./src/python/scripts/resize_ss.py`

- **Step** : prealign slide seq to nissl using slicer 3D+IGT and export transform slide_id_t3a in txt format

- **Step** : create transformed image using script `./src/python/scripts/apply_t3a.py`

- **Step** : register prealigned slide seq image to nissl using slicer 3D+IGT and save fiducials slide_id_f and slide_id_t (save as mark ups fiducial and change extension to csv)

- **Step** : run python script to get list of region ids - provide path to transforms and input, points, and nrrd file, and slide_index for transform, slice index in nrrd file, output list of region ids
