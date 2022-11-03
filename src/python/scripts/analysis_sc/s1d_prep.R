#ifndef S1D_PREP_R
#define S1D_PREP_R

# R script used as subroutine to perform regionwise aggregation
#  
# Created by Mukund on 2022-10-28


library(anndata)
Sys.setenv(RETICULATE_PYTHON = "/Users/mraj/opt/miniconda3/envs/Rv4/bin/python3")
library(reticulate)
library(Matrix)
library(Matrix.utils)

args = commandArgs(trailingOnly=TRUE)
print(args)

ip_folder_labels <- args[1]
ip_folder_counts <- args[2]
op_folder <- args[3]
# start_pid <- as.integer(args[4])
# end_pid <- as.integer(args[5])
pid <- as.integer(args[4])

print(paste('R processing pid:', pid))
# pid_idx <- pid_idx + 1
apid = pid
if (pid==5 || pid==77 || pid==167)
  apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

# get data matrix of puck
ip_counts_file = paste0(ip_folder_counts, "/dd", apid, "_CTMapping.h5ad")
# print(ip_counts_file)
ad <- read_h5ad(ip_counts_file)

# get labels factor
ip_labels_file = paste0(ip_folder_labels, '/agg_labels_',sprintf("%03d", apid),'.csv')
# print(ip_labels_file)
labels <- factor(read.csv(file=ip_labels_file, header=FALSE)$V1)

# print(length(labels))

# filter labels to keep only labels for selected beads
selected_ids <- ad$var$rowName
# browser()
# print(c('selected_ids length', length(selected_ids)))
# print(c(min(selected_ids), max(selected_ids)))

selected_ids <- selected_ids + 1 # adjust for R indexing
selectedBeadlabels <- labels[selected_ids]



# transpose X

X <- t(ad$X)
# X <- ad$X

# print(length(selectedBeadlabels))
# print(dim(X))
# print(c('min selected ids', min(selected_ids), max(selected_ids)))
# print(tail(ad$var,10))

## temp start

# ip_pos_file = paste0('/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds/', '/allen_anno_data_',sprintf("%03d", apid),'.csv')
# print(ip_labels_file)
# positions <- read.csv(file=ip_pos_file, header=FALSE, stringsAsFactors = TRUE)
# print(tail(positions[selected_ids,], 10))
# print(tail(selected_ids, 10))
# print ('end')
## temp ends

b <- aggregate.Matrix(X, selectedBeadlabels, FUN=sum()) # regions x chosen_cells(s)

print(dim(b))

ad <- AnnData(
  X = b
  #obs = data.frame(group = c("a", "b"), row.names = c("s1", "s2")),
  #var = data.frame(type = c(1L, 2L), row.names = c("var1", "var2")),
  
)

op_file <- paste0(op_folder, '/aggr_scores_', sprintf("%03d", apid), ".h5ad")
write_h5ad(ad, op_file);

#bead_counts = table(labels)
#op_file2 <-paste0(io_folder_interim, '/aggr_num_beads_', sprintf("%03d", apid), ".csv")
#write.table(bead_counts, file=op_file2,sep=',', quote=FALSE, row.names = FALSE)

gc()

#endif /* S1D_PREP_R */
