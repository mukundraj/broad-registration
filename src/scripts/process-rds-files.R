## Main file for processing slide-seq RDS files in a folder and creating 
## PNG images
## Created on 01-Dec-2021 by Mukund Raj



source('src/utils/jonah_utils.R')

# Input folder path

# Output folder path



mega_puck = mcreadRDS("../../../data/MBASS_1_For_Mukund/02_RDSs/Puck_210203_12_RCTD.RDS")

object_seurat = mcreadRDS("../../../data/MBASS_1_For_Mukund/02_RDSs/Puck_210203_12_Seurat.RDS")

# Read list of file names


# Process each file and save output



