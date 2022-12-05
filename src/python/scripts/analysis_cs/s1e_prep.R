# Prepares regionwise aggregated, greater than threshold bead scores for each
# mapped celltype. Similar to ../analysis/s9d_prep.R but with an additional
# clamping step. Also, similar to analysis/s9g_nz_scores.py for genes.
#
# Called as subroutine by analysis_cs/s2c_nz_scores.py
#
# Created by Mukund on 2022-12-04 

library(anndata)
Sys.setenv(RETICULATE_PYTHON = "/Users/mraj/opt/miniconda3/envs/Rv4/bin/python3")
library(reticulate)
library(Matrix)
library(Matrix.utils)

args = commandArgs(trailingOnly=TRUE)
print(args)

ip_folder_labels <- args[1]
ip_folder_scores <- args[2]
io_folder_interim <- args[3]
pid <- as.integer(args[4])

print(paste('R processing pid:', pid))
apid = pid
if (pid==5 || pid==77 || pid==167)
  apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

# get data matrix of puck
ip_scores_file = paste0(ip_folder_scores, "/dd", apid, "_CTMapping.h5ad")
# print(ip_scores_file)
ad <- read_h5ad(ip_scores_file)

# get labels factor
ip_labels_file = paste0(ip_folder_labels, '/agg_labels_',sprintf("%03d", apid),'.csv')
# print(ip_labels_file)
labels <- factor(read.csv(file=ip_labels_file, header=FALSE)$V1)

selected_ids <- ad$var$rowName
selected_ids <- selected_ids + 1 # adjust for R indexing
selectedBeadlabels <- labels[selected_ids]

X <- t(ad$X)

# making data binary by clamping upper bound of expression count to 1
threshold <- 0.3

ind <- X@x >= threshold
X@x[ind] <- 1 # set to 1 if greater than threshold to consider bead as containing celltype

X@x[X@x<threshold]<-0 # set all values less than threshold 0.3 to 0

print (dim(X))
print(length(selectedBeadlabels))
# aggregating regionwise
b <- aggregate.Matrix(X,selectedBeadlabels, FUN=sum()) # regions x mapped_celltypes(s)


ad <- AnnData(
  X = b
  #obs = data.frame(group = c("a", "b"), row.names = c("s1", "s2")),
  #var = data.frame(type = c(1L, 2L), row.names = c("var1", "var2")),
  
)

op_file <- paste0(io_folder_interim, '/nz_aggr_scores_', sprintf("%03d", apid), ".h5ad")
write_h5ad(ad, op_file);
print('wrote')

# bead_scores = table(labels)
bead_scores = table(selectedBeadlabels)
op_file2 <-paste0(io_folder_interim, '/nz_aggr_num_beads_', sprintf("%03d", apid), ".csv")
write.table(bead_scores, file=op_file2,sep=',', quote=FALSE, row.names = FALSE)

gc()

# https://stackoverflow.com/questions/33775291/r-matrix-set-particular-elements-of-sparse-matrix-to-zero
