
library(anndata)
Sys.setenv(RETICULATE_PYTHON = "/Users/mraj/opt/miniconda3/envs/Rv4/bin/python3")
library(reticulate)
library(Matrix)
library(Matrix.utils)


ip_folder_labels <- '/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s9_analysis/s9d/region_labels_cshl'
ip_folder_scores <- '/Users/mraj/Desktop/work/data/mouse_atlas/cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW'
# io_folder_interim <- args[3]
pid <- 41

apid <- pid

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
threshold <-0.30

ind <- X@x >= threshold
X@x[ind] <- 1 # set to 1 if greater than threshold to consider bead as containing celltype
X@x[X@x<threshold]<-0 # set all values less than threshold 0.3 to 0

# print (dim(X))kk
# save(X, selectedBeadlabels, file='./X.RData')
# print(length(selectedBeadlabels))
# aggregating regionwise
b <- aggregate.Matrix(X,selectedBeadlabels, FUN=sum()) # regions x mapped_celltypes(s)
dim(b)

which(b['961',]>0)
write.csv(X[,1927], "temp.csv")

write.csv(selectedBeadlabels, "temp1.csv")




