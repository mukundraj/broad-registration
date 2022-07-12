
library(anndata)
Sys.setenv(RETICULATE_PYTHON = "/Users/mraj/opt/miniconda3/envs/Rv4/bin/python3")
library(reticulate)
library(Matrix)
library(Matrix.utils)

args = commandArgs(trailingOnly=TRUE)
print(args)

#ip_file <- '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s8_raw_data/integrated_mats/ad_counts_21.h5ad'

#ad <- read_h5ad(ip_file)

# ip_folder_labels <- "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/region_labels"
# ip_folder_counts <- "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s8_raw_data/integrated_mats"
# io_folder_interim <- "/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/interim"
# start_pid <- 1
# end_pid <- 2

ip_folder_labels <- args[1]
ip_folder_counts <- args[2]
io_folder_interim <- args[3]
# start_pid <- as.integer(args[4])
# end_pid <- as.integer(args[5])
pid <- as.integer(args[4])

# pucks <- seq(start_pid,end_pid,2)

# pid_idx <- 0
# for (pid in pucks) {
  print(paste('R processing pid:', pid))
  # pid_idx <- pid_idx + 1
  apid = pid
  if (pid==5 || pid==77 || pid==167)
    apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment
  
  # get data matrix of puck
  ip_counts_file = paste0(ip_folder_counts, "/ad_counts_", apid, ".h5ad")
  # print(ip_counts_file)
  ad <- read_h5ad(ip_counts_file)
  
  # get labels factor
  ip_labels_file = paste0(ip_folder_labels, '/agg_labels_',sprintf("%03d", apid),'.csv')
  # print(ip_labels_file)
  labels <- factor(read.csv(file=ip_labels_file, header=FALSE)$V1)
  
  X <- t(ad$X)
  b <- aggregate.Matrix(X,labels, FUN=sum()) # regions x chosen_gene(s)
  
  
  ad <- AnnData(
    X = b
    #obs = data.frame(group = c("a", "b"), row.names = c("s1", "s2")),
    #var = data.frame(type = c(1L, 2L), row.names = c("var1", "var2")),
    
  )
  
  op_file <- paste0(io_folder_interim, '/aggr_counts_', sprintf("%03d", apid), ".h5ad")
  write_h5ad(ad, op_file);
  
  bead_counts = table(labels)
  op_file2 <-paste0(io_folder_interim, '/aggr_num_beads_', sprintf("%03d", apid), ".csv")
  write.table(bead_counts, file=op_file2,sep=',', quote=FALSE, row.names = FALSE)
  
  gc()
# }
