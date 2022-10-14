"""
Read and process bead x cell score matrices.


Usage:

Usage example:

python src/python/scripts/analysis_sc/s1b_beadxcell.py

Supplementary:

gsutil -m rsync -r gs://macosko_data/jlanglie/scp/03_All_MBASS_Mapping_Mega_Matrix ./03_All_MBASS_Mapping_Mega_Matrix

Created by Mukund on 2022-10-13

"""

from produtils import dprint
import anndata as ann

in_folder = "/Users/mraj/Desktop/work/data/mouse_atlas/single_cell/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix/"

inp_file = f'{in_folder}/dd1_CTMapping.h5ad'

data = ann.read_h5ad(inp_file)

dprint(data)
dprint(data.X)
