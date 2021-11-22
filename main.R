library(stringr)
library(magrittr)
library(dplyr)
library(pals)
library(ggplot2)
library(glue)
library(purrr)



# Plot GAPDH for mukund

# Helps remove outliers
clamp <- function(lst, max = Inf, min = -max, rm = FALSE) {
  if (rm == TRUE) {
    lst = lst[lst < max]
    lst = lst[lst > min]
    return(lst)
  }
  else {
    return(pmax(pmin(lst, max), min))
  }
}

# Faster import/export
mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}
mcreadRDS <- function(file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -d -c -p",mc.cores," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

# Import file
mega_puck = mcreadRDS("/home/mraj/data/forMukund/Puck_210203_04.RDS")

# If want to plot GAPDH, find the gene in the gene names
grep("gapdh", rownames(mega_puck), value=T, ignore.case = T)
# Gapdh

# Optional but remove outright beads with less than 150 UMU's
UMI_cutoff_bcds = names(which(
  colSums(mega_puck@assays$Spatial@counts) > 150
))

# # (1) Plot GAPDH (or any gene)
# gene_val = mega_puck@assays$Spatial@counts["Gapdh", UMI_cutoff_bcds]
# # If want log. Note that log(0) is bad, so replace those with 0
# # gene_val = log10(gene_val) %>% {.[is.infinite(.)]=0; .}

# (2) Plot UMI counts per bead
gene_val = colSums(mega_puck@assays$Spatial@counts[, UMI_cutoff_bcds])
# Remove outliers
gene_val_clamp =  gene_val %>% clamp(max=8000)
hist(gene_val_clamp)
# Optionally log10

# # (3) Plot number of non-zero genes per bead
# gene_val = colSums(mega_puck@assays$Spatial@counts[, UMI_cutoff_bcds] > 0)
# hist(gene_val)
# # Outliers
# gene_val_clamp =  gene_val %>% clamp(max=5000)
# hist(gene_val_clamp)
# # Optional log10

coords = mega_puck@images$image@coordinates[UMI_cutoff_bcds, ]

load("tmp_data/coords.Rdata")
load("tmp_data/gene_val.Rdata")
source("plot_puck.R")
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
ggsave("output/ss.png", plot=last_plot(), dpi=96, scale=8, limitsize=FALSE)
