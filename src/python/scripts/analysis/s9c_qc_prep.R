#!/usr/bin/env Rscript
#
# Reads gene and output nRegion x nPuck matrix for genes specified in argument
#
# Usage example:
#
#  Rscript src/python/scripts/analysis/s9c_qc_prep.R
#
# Created by Mukund on 2022-06-11

library(anndata)
Sys.setenv(RETICULATE_PYTHON = "/Users/mraj/opt/miniconda3/envs/Rv4/bin/python3")
library(reticulate)
library(Matrix)
library(Matrix.utils)



args = commandArgs(trailingOnly=TRUE)
print(args)
# 
# gene = 'Ndnf'
# gene = 'Gad1'
# ip_folder_labels <- "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/aggregated_labels"
# ip_folder_counts <- "/Users/mraj/Desktop/work/data/temp_data/2022-05-04/integrated_mats"
# r_op_filename <- "/Users/mraj/Desktop/work/data/mouse_atlas/temp/Ndnf.csv"
# region_ids <- '1089'
# start_pid <- 79
# end_pid <- 80
# nonzero <- FALSE
# region_ids <- '1;2;3;4;5'
# region_ids <- '1006;670;1086;1111;9;461;1089'

gene <- args[1]
ip_folder_labels <- args[2]
ip_folder_counts <- args[3]
r_op_filename <- args[4]
region_ids <- args[5]
start_pid <- as.integer(args[6])
end_pid <- as.integer(args[7])
nonzero <- args[8]=='True'
print(start_pid)
print(end_pid)
print(region_ids)
print(nonzero)

br_layer_ids <- strsplit(region_ids,";")

rids_lists <- vector(mode = "list", length = length(br_layer_ids[[1]]))
ii=0
for (brlstr in br_layer_ids[[1]]){
  ii=ii+1
  
  brl_ids <- strsplit(brlstr[[1]], ',')
  rids_lists[[ii]] <- vector(mode = "list", length = length(brl_ids[[1]]))
  # brl_ids <- as.numeric(unlist(brl_ids))
  # brl_ids <- unlist(brl_ids)

  jj=0
  for (brl_id in brl_ids){
        jj=jj+1
        
        rids_lists[[ii]][[jj]] <- brl_id
  }
}

pucks <- seq(start_pid,end_pid,2)
gene_mat <- matrix(0, nrow = length(rids_lists), ncol = length(pucks))
print(pucks)
# loop over pucks
pid_idx <- 0
for (pid in pucks) {
  print(pid)
  pid_idx <- pid_idx + 1
  apid = pid
  if (pid==5 || pid==77 || pid==167)
    apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment
   


  # get data matrix of puck
  ip_counts_file = paste0(ip_folder_counts, "/ad_counts_", apid, ".h5ad")
  print(ip_counts_file)
  ad <- read_h5ad(ip_counts_file)
  X <- t(ad$X)
  print(dim(X))
  print(length(ad$obs_names))
  # filter selected gene
  gene_idx <- ad$obs_names==gene
  X_filt <- X[,gene_idx]
  totals_umis <- rowSums(X) # total umis f
  total_umis_gene <- totals_umis[gene_idx]
  #print(length(totals_umis))
  print(c('len X_filt', length(X_filt)))
  
  
  
  # get labels factor
  ip_labels_file = paste0(ip_folder_labels, '/agg_labels_',sprintf("%03d", apid),'.csv')
  print(ip_labels_file)
  labels <- factor(read.csv(file=ip_labels_file, header=FALSE)$V1)
  
  # aggregate over regions
  b<-aggregate.Matrix(X_filt,labels, FUN=sum())
  
  cmat<-aggregate.Matrix(X,labels, FUN=sum())
  #print(b)
  print('entering loop')
  # populate heatmap column for this puck by iterating over required regions
  for (ii in seq(length(rids_lists))){ # looping over region
   
    sum <- 0;
    for (jj in seq(length(rids_lists[[ii]]))){ ## looping over sub region
      
      # print("jj")
      # print(rids_lists[ii][jj][[1]][[1]]);
      # print(levels(labels))
      # print(any(levels(labels)==rids_lists[ii][jj][[1]][[1]]))
      # print(any(levels(labels)=="1089"))
      
      
      # print('total_region_umis')
      # print(total_region_umis)
      
      # get gene count in sub region
      if (any(levels(labels)==rids_lists[ii][jj][[1]][[1]])){
        total_region_umis <- sum(cmat[rids_lists[ii][jj][[1]][[1]],])
        cnt <- b[rids_lists[ii][jj][[1]][[1]], , drop = FALSE]@x
        if (length(cnt)==1){
          # update count in mat
          if (nonzero){
            gene_mat[ii, pid_idx] <- gene_mat[ii, pid_idx]+1
          }else{
            gene_mat[ii, pid_idx] <- gene_mat[ii, pid_idx]+cnt/total_region_umis
            #print(c(ii,jj,rids_lists[[ii]][[jj]][[1]][[1]]))
          }
          
        }
       
      }
    }
  }
  
  gc()
  
}

# write out output matrix csv
write.table(gene_mat, file=r_op_filename, row.names=FALSE, col.name=FALSE, sep=',')

# matrix(, nrow = 2, ncol = 2)
# any(levels(labels)=="1016")
# b['1016', , drop = FALSE]@x
# write.table(mm, file='temp.csv', row.names=FALSE, col.name=FALSE, sep=',')
# 
# https://stat.ethz.ch/pipermail/r-help/2011-April/275773.html
# r_op_filename = '/Users/mraj/Desktop/work/data/mouse_atlas/temp/Ndnf-cnt.csv'
# gene_mat <- as.matrix(read.table(file=r_op_filename, sep=','))
# colnames(gene_mat) = pucks
# rownames(gene_mat) = c('L1', 'L2/3','L4', 'L5', 'L6', 'L7', 'Hippo Formation' )
# heatmap(gene_mat, Colv = NA, Rowv = NA, scale="column")

# #Initial Heatmap 
# install.packages("tidyverse")
# library(tidyverse)
# gene_mat %>% 
#   as.data.frame() %>%
#   rownames_to_column("f_id") %>%
#   pivot_longer(-c(f_id), names_to = "samples", values_to = "counts") %>%
#   mutate(samples= fct_relevel(samples,colnames(dat))) %>%
#   ggplot(aes(x=samples, y=f_id, fill=counts)) + 
#   geom_raster() + 
#   scale_fill_viridis_c()
