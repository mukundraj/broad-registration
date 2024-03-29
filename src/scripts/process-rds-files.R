## Main file for processing slide-seq RDS files in a folder and creating 
## PNG images
## Created on 01-Dec-2021 by Mukund Raj

source('src/utils/jonah_utils.R')
source('src/utils/plot_puck.R')
source('src/utils/mapper.R')

Sys.setenv(R_CONFIG_ACTIVE = "current")

data_path <- config::get("mapping_data") # data sheet mapping slide seq to nissl

# data_path <- 'input/mapping_data.tsv' # data sheet mapping slide seq to nissl
mapper <- map_slseq_nissl(data_path)


# Input folder path
# input_path = "../../../data/MBASS_1_For_Mukund/02_RDSs/sdata"
# input_path <- "/Users/mraj/Desktop/work/data/forMukund/allRDSs"

input_path <- config::get("input_path")

files = list.files(input_path, full.names=TRUE)

# Output folder path
#op_path = "output/ss/nUMI/"
# op_path <- "output/ss/nUMI_v2/"
op_path <- config::get("op_path")


## Sample code to read rctd and seurat objects
# mega_puck = mcreadRDS("../../../data/MBASS_1_For_Mukund/02_RDSs/Puck_210203_12_RCTD.RDS")
# object_seurat = mcreadRDS("../../../data/MBASS_1_For_Mukund/02_RDSs/Puck_210203_12_Seurat.RDS")

# Read list of file names

for (i in 1:length(files))
# for (i in 1:5)
{
  nissl_id <- mapper[basename(files[i])]
  #nissl_id <- nissl_id[[1]]*2 - 1 # to align with Chuck's numbering scheme
  # if (nissl_id!=219) # sample dataset id
  #   next
  
 
  
  print (paste(i, nissl_id, "Processing", files[i]))
  # if (nissl_id==143){
  #   print (paste(i, nissl_id, "Processing", files[i]))
  # }
  # else{
  #   next
  # }
  
  
 
  object_seurat = mcreadRDS(files[i])
  
  # Optional but remove outright beads with less than 150 UMU's
  UMI_cutoff_bcds = names(which(
    # colSums(object_seurat@assays$Spatial@counts) > 150
    colSums(object_seurat@assays$Spatial@counts) > -1
  ))
  
  coords = object_seurat@images$image@coordinates[UMI_cutoff_bcds, ]
  write.table(coords,paste("output/ss/bead_coords/coords_", nissl_id, ".csv", sep=""), sep=",", row.names = FALSE, col.names=FALSE)
  
  
  gene_val = colSums(object_seurat@assays$Spatial@counts[, UMI_cutoff_bcds])
  gene_val_clamp =  gene_val %>% clamp(max=8000)
  
# Note clamp doesn't remove, just clamps. Playing around fast, I like 0.95 - Jonah
  clamp(gene_val, quantile(gene_val, 0.999)) %>% {plot(density(.))}
  gene_val_clamp =  gene_val %>% {clamp(., max=quantile(., 0.99))}
  gene_val_clamp  = log10(gene_val)
  gene_val_clamp  = log10(gene_val_clamp)
  
 
  # new extent code starts
  
  minmax_x = c(min(coords$x), max(coords$x))
  minmax_y = c(min(coords$y), max(coords$y))
  
  print(paste("minmax_x", minmax_x))
  print(paste("minmax_y", minmax_y))
  xvals <- c(minmax_x[1], minmax_x[2], minmax_x[1], minmax_x[2])
  yvals <- c(minmax_y[2], minmax_y[1], minmax_y[1], minmax_y[2])

  corners <- data.frame(x=xvals, y=yvals)
 
  print(corners)
  
  # new extent code ends
  
  
  gene_val_cat = rep_along(gene_val_clamp, "FALSE")
  gene_val_cat = c(gene_val_cat, "TRUE", "TRUE", "TRUE", "TRUE")
  gene_val_cat  = factor(gene_val_cat)
  plot_puck(val = gene_val_cat, #gene_val_clamp,
            # corners = corners,
            corners = NULL,
            # Could also pass in puck = mega_puck but wanted to subset UMI_cutoff
            coords = rbind(coords, corners),
            # Don't make vector with thousands of circles
            raster=F,
            # Custom color scheme with alpha. Going from transparent blue -> opaque yellow
            # Orders by valye ut still nice to push down the opacity of background beads
            # pal=colorRampPalette(c(alpha("blue", 0.01),
            #                        alpha("green", 0.2),
            #                        alpha("yellow", 0.7)), alpha=T)(100)
            # pal=colorRampPalette(c(alpha("white", 0.01),
            #                        alpha("black", 1)), alpha=T)(100)
            pal = c("#AAAAAA", "RED")
  )




  print(paste('output:', op_path, nissl_id))
  ggsave(paste(op_path,nissl_id,".png", sep=""), plot=last_plot(), dpi=72, scale=8, limitsize=FALSE, height = 10, width = 10)
  
}

## command stash
# fd seurat -x cp {} sdata/
# fd rctd -x cp {} rdata/
  
