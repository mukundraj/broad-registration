## Main file for processing slide-seq RDS files in a folder and creating 
## PNG images
## Created on 01-Dec-2021 by Mukund Raj

source('src/utils/jonah_utils.R')
source('src/utils/plot_puck.R')
source('src/utils/mapper.R')

Sys.setenv(R_CONFIG_ACTIVE = "current")

data_path <- 'input/mapping_data.tsv' # data sheet mapping slide seq to nissl
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

# for (i in 1:length(files))
for (i in 1:2)
{
  nissl_id <- mapper[basename(files[i])]
  nissl_id <- nissl_id[[1]]*2 - 1 # to align with Chuck's numbering scheme
  # if (nissl_id!=219) # sample dataset id
  #   next

  
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
  
  
  
  
  
  ggsave(paste(op_path,nissl_id,".png", sep=""), plot=last_plot(), dpi=96, scale=8, limitsize=FALSE, height = 10, width = 10)
  
}


## command stash
# fd seurat -x cp {} sdata/
# fd rctd -x cp {} rdata/
  
