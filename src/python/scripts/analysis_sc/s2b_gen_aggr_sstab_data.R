# R script to generate zarr file for SingleCell tab data aggregated by clades and cell classes
#
#
# Supplementary:
#
# gsutil rsync -r ~/Desktop/work/data/mouse_atlas/single_cell/s2/agged_zarr gs://bcdportaldata/batch_231112/single_cell
#
# Created by Mukund on 2023-11-13
#
# Params:
#
# avg_loc - location of avg expression values 
# nzpct_loc - location of nz_pct (percent nonzero) values
# cladeloc - location of clades info

avg_loc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_321017.zarr/avg/X'
nzpct_loc  <-  '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_321017.zarr/nz_pct/X'
cladeloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_321017.zarr/metadata/clades'
cellclassloc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_321017.zarr/metadata/cellclasses'
op_loc  <- '~/Desktop/work/data/mouse_atlas/single_cell/s2/agged_zarr/'

# Params end



# zloc <- '/Users/mraj/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_321017.zarr'
# zloc <- '/Users/mraj/Desktop/work/data/mouse_atlas/single_cell/s1/scZarr_321017.zarr/counts/X'


library(Rarr)
library(Matrix)
library(Matrix.utils)


# read_array_metadata(zloc)

# avg_array  <- zarr_overview(zloc)
avg_array  <- read_zarr_array(avg_loc)

clades  <- read_zarr_array(cladeloc)
# convert clades to list
clades_factor  <- factor(clades)


# aggregate by zo by clades

agg_avg  <- aggregate.Matrix(avg_array, clades_factor, FUN = 'mean')
agg_avg_array  <- as.matrix(agg_avg)

# # write to file
chunk_dim = c(dim(agg_avg)[1], 1)
create_empty_zarr_array(op_loc, data_type="double", dim=dim(agg_avg), chunk_dim=chunk_dim)

agged_avgs_loc  <- file.path(op_loc, "aggedSCdata.zarr/byclades/avg/X")
# create folder if not exists
dir.create(agged_avgs_loc, recursive = TRUE, showWarnings = FALSE)
write_zarr_array( x=agg_avg_array, zarr_array_path=agged_avgs_loc, chunk_dim=chunk_dim,)

# aggregage nz_pct by clades
agg_nz_pct  <- aggregate.Matrix(read_zarr_array(nzpct_loc), clades_factor, FUN = 'mean')
agg_nz_pct_array  <- as.matrix(agg_nz_pct)

# write to file
chunk_dim = c(dim(agg_nz_pct)[1], 1)
agged_nz_pct_loc  <- file.path(op_loc, "aggedSCdata.zarr/byclades/nz_pct/X")
# create folder if not exists
dir.create(agged_nz_pct_loc, recursive = TRUE, showWarnings = FALSE)
write_zarr_array( x=agg_nz_pct_array, zarr_array_path=agged_nz_pct_loc, chunk_dim=chunk_dim,)


## aggregate avg by cell classes

cellclasses  <- read_zarr_array(cellclassloc)
# convert clades to list
cellclasses_factor  <- factor(cellclasses)

agg_avg  <- aggregate.Matrix(avg_array, cellclasses_factor, FUN = 'mean')
agg_avg_array  <- as.matrix(agg_avg)
  
# write to file
chunk_dim = c(dim(agg_avg)[1], 1)
agged_avgs_loc  <- file.path(op_loc, "aggedSCdata.zarr/bycellclasss/avg/X")
# create folder if not exists
dir.create(agged_avgs_loc, recursive = TRUE, showWarnings = FALSE)
write_zarr_array( x=agg_avg_array, zarr_array_path=agged_avgs_loc, chunk_dim=chunk_dim,)


# aggregage nz_pct by cell classes
agg_nz_pct  <- aggregate.Matrix(read_zarr_array(nzpct_loc), cellclasses_factor, FUN = 'mean')
agg_nz_pct_array  <- as.matrix(agg_nz_pct)

# write to file
chunk_dim = c(dim(agg_nz_pct)[1], 1)
agged_nz_pct_loc  <- file.path(op_loc, "aggedSCdata.zarr/bycellclasss/nz_pct/X")
# create folder if not exists
dir.create(agged_nz_pct_loc, recursive = TRUE, showWarnings = FALSE)
write_zarr_array( x=agg_nz_pct_array, zarr_array_path=agged_nz_pct_loc, chunk_dim=chunk_dim,)


print("Done")
