""" 

python src/python/scripts/analysis_cs/debug.py \

"""

import anndata as ann
from scipy.sparse import csr_matrix
import numpy as np
from produtils import dprint
import csv

d = ann.read_h5ad('/Users/mraj/Desktop/work/data/mouse_atlas/cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW/dd41_CTMapping.h5ad')

# tmp = colnames(b)
# which(tmp=='301-1-0=Micro_Selplg_Siglech')

# cellname = '301-1-0=Micro_Selplg_Siglech'
cellname = '300-0-0-1-0-0-1-0-2-1-1=Ex_Rorb_Ptpn20'
tmp = d.obs.axes[0].tolist()
idx = tmp.index(cellname)
dprint(cellname, idx)

d.obs.iloc[[idx]] # 1926 is Micro_Selplg_Siglech
row = csr_matrix(d.X.getrow(idx)).todense()
ninds = np.where(row>=0.3)[1]

# inds in all beads file for this celltype
oinds = [ (int(d.var.rowName.iloc[i]),i) for i in ninds] # orig ind, jonah anndata ind

label_data_folder = f'/Users/mraj/Desktop/work/data/mouse_atlas/v3/s2/bead_ccf_labels_allbds'
# inds in periform area
nis_id_str = str(41).zfill(3)
labels_csv_file = f'{label_data_folder}/allen_anno_data_{nis_id_str}.csv'
# dprint(labels_csv_file)
region_names = []
with open(labels_csv_file, newline='\n') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        region_names.append(row[4])


finds = [ (i,i2) if region_names[i]=='Piriform area' else 'XXX' for (i,i2) in oinds ] # filtered inds

out_count = finds.count('XXX')
tot_count = len(finds)
in_frac = (tot_count-out_count)/tot_count

dprint('out_count', out_count)
dprint('tot_count', tot_count)
dprint('in_frac', in_frac)



# dprint(finds)


# for i in finds:
#     if i!='':
#         dprint(i, d.X[idx, i[1]])



