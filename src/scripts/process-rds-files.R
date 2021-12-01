## Main file for processing slide-seq RDS files in a folder and creating 
## PNG images
## Created on 01-Dec-2021 by Mukund Raj



source('src/utils/jonah_utils.R')



# Input folder path
input_path = "../../../data/MBASS_1_For_Mukund/02_RDSs/sdata"
files = list.files(input_path, full.names=TRUE)

# Output folder path


## Sample code to read rctd and seurat objects
# mega_puck = mcreadRDS("../../../data/MBASS_1_For_Mukund/02_RDSs/Puck_210203_12_RCTD.RDS")
# object_seurat = mcreadRDS("../../../data/MBASS_1_For_Mukund/02_RDSs/Puck_210203_12_Seurat.RDS")

# Read list of file names

# for (i in 1:length(files))
for (i in 1:4)
{
  print (paste("Processing", files[i]))
  object_seurat = mcreadRDS(files[i])
  
  # Optional but remove outright beads with less than 150 UMU's
  UMI_cutoff_bcds = names(which(
    colSums(object_seurat@assays$Spatial@counts) > 150
  ))
  
  coords = object_seurat@images$image@coordinates[UMI_cutoff_bcds, ]
  
  
  gene_val = colSums(object_seurat@assays$Spatial@counts[, UMI_cutoff_bcds])
  gene_val_clamp =  gene_val %>% clamp(max=8000)
  
  
  plot_puck(val = gene_val_clamp,
            # Could also pass in puck = mega_puck but wanted to subset UMI_cutoff
            coords = coords,
            # Don't make vector with thousands of circles
            raster=F,
            # Custom color scheme with alpha. Going from transparent blue -> opaque yellow
            # Orders by valye ut still nice to push down the opacity of background beads
            pal=colorRampPalette(c(alpha("blue", 0.01),
                                   alpha("green", 0.2),
                                   alpha("yellow", 0.7)), alpha=T)(100)
  )
  
  ggsave(paste("output/ss/nUMI/",i,".png", sep=""), plot=last_plot(), dpi=96, scale=8, limitsize=FALSE)
  
}


# Process each file and save output





## command stash
# fd seurat -x cp {} sdata/
# fd rctd -x cp {} rdata/
  
