# Get clade/class matrices with both avg and pct exprs



# Created by Mukund on 2024-09-27

genes = read.csv('genes.csv', header=FALSE)
colnames(a) <-genes[,1]

##################################################
# directory mouse_atlas/single_cell/s2/agged_h5ad/
##################################################

avg_class <- read_h5ad('cellclasses_agged_avgs.h5ad')
avg_clade <- read_h5ad('clades_agged_avgs.h5ad')
pct_class <- read_h5ad('cellclasses_agged_nz_pct.h5ad')
pct_clade <- read_h5ad('clades_agged_nz_pct.h5ad')

genes = read.csv('genes.csv', header=FALSE)

colnames(avg_class) <-genes[,1]
colnames(avg_clade) <-genes[,1]
colnames(pct_class) <-genes[,1]
colnames(pct_clade) <-genes[,1]

qsave(avg_class$X, 'avg_class.qs')
qsave(avg_clade$X, 'avg_clade.qs')
qsave(pct_class$X, 'pct_class.qs')
qsave(pct_clade$X, 'pct_clade.qs')
