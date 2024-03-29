# R script to generate anndata files for SingleCell tab data aggregated by clades and cell classes - this data is reformatted to zarr format in s3a_gen_aggr_sstab_zarr.py. Also creates cladeAnnotations.csv for downstream processing by s3a_gen_aggr_sstab_zarr.py
#
# Example:
#
# Rscript ./src/python/scripts/analysis_sc/s2b_gen_aggr_sstab_h5ad.R
#
# Supplementary:
#
# gsutil rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s2/agged_zarr gs://bcdportaldata/batch_231112/single_cell
# gsutil rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s2/agged_zarr/aggedSCdata.zarr/bycellclasss/counts gs://bcdportaldata/batch_231112/single_cell/bycellclass/counts
#
# Created by Mukund on 2023-11-13
#
# Params:
#
# avg_loc - location of avg expression values 
# nzpct_loc - location of nz_pct (percent nonzero) values
# cladeloc - location of clades info

# avg_loc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_231207.zarr/avg/X'
# nzpct_loc  <-  '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_231207.zarr/nz_pct/X'
# counts_loc <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_231207.zarr/counts/X'
# cladeloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_231207.zarr/metadata/clades'
# cladeAnnotationsloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_231207.zarr/metadata/cladeAnnotations'
# cellclassloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_231207.zarr/metadata/cellclasses'

avg_loc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr/avg/X'
nzpct_loc  <-  '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr/nz_pct/X'
counts_loc <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr/counts/X'
cladeloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr/metadata/clades'
cladeAnnotationsloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr/metadata/cladeAnnotations'
cellclassloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_240318.zarr/metadata/cellclasses'
op_loc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s2/agged_h5ad'
op_cladeanno_loc <- '~/Desktop/work/data/mouse_atlas/single_cell/s2/cladeAnnotations.csv'

# Params end


library(Rarr)
library(Matrix)
library(Matrix.utils)
library(anndata)


# read_array_metadata(zloc)

clades  <- read_zarr_array(cladeloc)
# convert clades to list
clades_factor  <- factor(clades)

cladeAnnotations  <- read_zarr_array(cladeAnnotationsloc)

# create a new dataframe with clades and cladeAnnotations as columns
cladeAnnotations  <- data.frame(clades, cladeAnnotations)

# write to csv
write.csv(cladeAnnotations, file = op_cladeanno_loc)

# avg_array  <- zarr_overview(zloc)
# avg_array  <- read_zarr_array(avg_loc)

## agg by clades
# aggregate avg by clades

agg_avg  <- aggregate.Matrix(read_zarr_array(avg_loc), clades_factor, fun = 'mean')
agg_avg_array  <- as.matrix(agg_avg)
ad <- AnnData(X = agg_avg_array, obs = data.frame(group = rownames(agg_avg_array), row.names = rownames(agg_avg_array)) )
agged_avgs_loc  <- file.path(op_loc, "/clades_agged_avgs.h5ad")
write_h5ad(ad, agged_avgs_loc, )


# # # write to file
# chunk_dim = c(dim(agg_avg)[1], 1)
# create_empty_zarr_array(op_loc, data_type="double", dim=dim(agg_avg), chunk_dim=chunk_dim)

# agged_avgs_loc  <- file.path(op_loc, "aggedSCdata.zarr/byclades/avg/X")
# # create folder if not exists
# dir.create(agged_avgs_loc, recursive = TRUE, showWarnings = FALSE)
# write_zarr_array( x=agg_avg_array, zarr_array_path=agged_avgs_loc, chunk_dim=chunk_dim,)

# aggregage nz_pct by clades
agg_nz_pct  <- aggregate.Matrix(read_zarr_array(nzpct_loc), clades_factor, fun = 'mean')
agg_nz_pct_array  <- as.matrix(agg_nz_pct)
ad <- AnnData(X = agg_nz_pct_array, obs = data.frame(group = rownames(agg_nz_pct_array), row.names = rownames(agg_nz_pct_array)) )
agged_nz_pct_loc  <- file.path(op_loc, "/clades_agged_nz_pct.h5ad")
write_h5ad(ad, agged_nz_pct_loc, )

# # write to file
# chunk_dim = c(dim(agg_nz_pct)[1], 1)
# agged_nz_pct_loc  <- file.path(op_loc, "aggedSCdata.zarr/byclades/nz_pct/X")
# # create folder if not exists
# dir.create(agged_nz_pct_loc, recursive = TRUE, showWarnings = FALSE)
# write_zarr_array( x=agg_nz_pct_array, zarr_array_path=agged_nz_pct_loc, chunk_dim=chunk_dim,)

## aggregate counts by clades

agg_counts = aggregate.Matrix(read_zarr_array(counts_loc), clades_factor, fun = 'sum')
agg_counts_array  <- as.matrix(agg_counts)
ad <- AnnData(X = agg_counts_array, obs = data.frame(group = rownames(agg_counts_array), row.names = rownames(agg_counts_array)) )
agged_counts_loc  <- file.path(op_loc, "/clades_agged_counts.h5ad")
write_h5ad(ad, agged_counts_loc, )

# # write to file
# chunk_dim = c(dim(agg_counts)[1], 1)
# agged_counts_loc  <- file.path(op_loc, "aggedSCdata.zarr/byclades/counts/X")
# # create folder if not exists
# dir.create(agged_counts_loc, recursive = TRUE, showWarnings = FALSE)
# write_zarr_array( x=agg_counts_array, zarr_array_path=agged_counts_loc, chunk_dim=chunk_dim,)

print("starting with cell classes")
# ## agg by cell classes
# # aggregate avg by cell classes

cellclasses  <- read_zarr_array(cellclassloc)
# convert clades to list
cellclasses_factor  <- factor(cellclasses)

agg_avg  <- aggregate.Matrix(read_zarr_array(avg_loc), cellclasses_factor, fun = 'mean')
agg_avg_array  <- as.matrix(agg_avg)
ad <- AnnData(X = agg_avg_array, obs = data.frame(group = rownames(agg_avg_array), row.names = rownames(agg_avg_array)) )
agged_avgs_loc  <- file.path(op_loc, "/cellclasses_agged_avgs.h5ad")
write_h5ad(ad, agged_avgs_loc, )
  
# # write to file
# chunk_dim = c(dim(agg_avg)[1], 1)
# agged_avgs_loc  <- file.path(op_loc, "aggedSCdata.zarr/bycellclasss/avg/X")
# # create folder if not exists
# dir.create(agged_avgs_loc, recursive = TRUE, showWarnings = FALSE)
# write_zarr_array( x=agg_avg_array, zarr_array_path=agged_avgs_loc, chunk_dim=chunk_dim,)


# aggregage nz_pct by cell classes
agg_nz_pct  <- aggregate.Matrix(read_zarr_array(nzpct_loc), cellclasses_factor, fun = 'mean')
agg_nz_pct_array  <- as.matrix(agg_nz_pct)
ad <- AnnData(X = agg_nz_pct_array, obs = data.frame(group = rownames(agg_nz_pct_array), row.names = rownames(agg_nz_pct_array)) )
agged_nz_pct_loc  <- file.path(op_loc, "/cellclasses_agged_nz_pct.h5ad")
write_h5ad(ad, agged_nz_pct_loc, )

# # write to file
# chunk_dim = c(dim(agg_nz_pct)[1], 1)
# agged_nz_pct_loc  <- file.path(op_loc, "aggedSCdata.zarr/bycellclasss/nz_pct/X")
# # create folder if not exists
# dir.create(agged_nz_pct_loc, recursive = TRUE, showWarnings = FALSE)
# write_zarr_array( x=agg_nz_pct_array, zarr_array_path=agged_nz_pct_loc, chunk_dim=chunk_dim,)

# # aggregate counts by cell classes

agg_counts = aggregate.Matrix(read_zarr_array(counts_loc), cellclasses_factor, fun = 'sum')
agg_counts_array  <- as.matrix(agg_counts)
ad <- AnnData(X = agg_counts_array, obs = data.frame(group = rownames(agg_counts_array), row.names = rownames(agg_counts_array)) )
agged_counts_loc  <- file.path(op_loc, "/cellclasses_agged_counts.h5ad")
write_h5ad(ad, agged_counts_loc, )

# # write to file
# chunk_dim = c(dim(agg_counts)[1], 1)
# agged_counts_loc  <- file.path(op_loc, "aggedSCdata.zarr/bycellclasss/counts/X")
# # create folder if not exists
# dir.create(agged_counts_loc, recursive = TRUE, showWarnings = FALSE)
# write_zarr_array( x=agg_counts_array, zarr_array_path=agged_counts_loc, chunk_dim=chunk_dim,)

print("Done")
