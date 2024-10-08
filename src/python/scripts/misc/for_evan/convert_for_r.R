# Converts ouptut of ./get_avg_mat.py to a qs file for use in R  
#
# Usage example: 
# 
# Rscript src/python/scripts/misc/for_evan/convert_for_r.R
#
# Created by Mukund on 2024-01-13

library(anndata)
library(qs)

ipfile <- "/Users/mraj/Desktop/work/data/temp1/avg_mat.h5ad"
ipfile_pct <- "/Users/mraj/Desktop/work/data/temp1/pct_mat.h5ad" 
ipfile_genes  <- "/Users/mraj/Desktop/work/data/temp1/genes.csv"
ipfile_cells  <- "/Users/mraj/Desktop/work/data/temp1/clusters.csv"

opfile <- "/Users/mraj/Desktop/work/data/temp1/avg_clusters.qs"
opfile_pct <- "/Users/mraj/Desktop/work/data/temp1/pct_clusters.qs"


# Read in the file
ip <- read_h5ad(ipfile)
ip_genes <- read.csv(ipfile_genes, header = FALSE)
ip_cells <- read.csv(ipfile_cells, header = FALSE)

# print(ip_genes)
# Convert to a dataframe with genes as columns and cells as rows
ip <- as.data.frame(as.matrix(ip), row.names = ip_cells$V1, col.names = ip_genes$V1)




colnames(ip) <- ip_genes$V1
# write to qs file
qsave(ip, opfile)


ip_pct <- read_h5ad(ipfile_pct)
ip_pct <- as.data.frame(as.matrix(ip_pct), row.names = ip_cells$V1, col.names = ip_genes$V1)
qsave(ip_pct, opfile_pct)

