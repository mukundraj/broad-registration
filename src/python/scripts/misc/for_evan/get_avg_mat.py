"""
This script reads the zarr file and writes the average matrix, genes and clusters to h5ad and csv files respectively.

Output directory is specified by op_dir in the beginning of the script.

Usage example:

python src/python/scripts/misc/for_evan/get_avg_mat.py \

Created by Mukund on 2024-01-13
"""


import zarr
import anndata as an
import os

ip_fname = "~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_231207.zarr"
op_dir = '~/Desktop/temp1'


# read zarr file
z = zarr.open(ip_fname, mode='r')


# create op_dir if not exists
if not os.path.exists(op_dir):
    os.makedirs(op_dir)


avg_mat = z.avg.X[:]

# convert to anndata
adata = an.AnnData(avg_mat)

# write to h5ad
adata.write(f'{op_dir}/avg_mat.h5ad')


genes = z.var.genes[:]

clusters = z.obs.clusters[:]

# keep only the part after '=' in cluster names
clusters = [cluster.split('=')[1] for cluster in clusters]

# write 
with open(f'{op_dir}/genes.csv', "w") as f:
    for gene in genes:
        f.write(gene+"\n")

with open(f'{op_dir}/clusters.csv', "w") as f:
    for cluster in clusters:
        f.write(cluster+"\n")




print("Done")


