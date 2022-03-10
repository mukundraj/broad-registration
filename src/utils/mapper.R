# Function to read in tsv data file with mapping between nissl and slide seq and 
# returning an associative array (List) with the mapping
# Created by Mukund Raj on 2022-02-02

# returns map between Jonah's Seurat RDS and nissl number
map_slseq_nissl <- function(data_path){

  # data <- read.table(file = data_path, sep = '\t', header = TRUE)
  data <- read.delim(data_path, sep="\t", header=TRUE, quote="\"")
  mapper = list()
  
  for (i in 1:length(data$Puck.ID)){
    key <- gsub(" ", "", paste(data$Puck.ID[i], "_Seurat.RDS", collapse=""))
    if(nchar(data$Puck.ID[i])>1){
      # print(data$Puck.ID[i])
      # print(data[[5]][i])
      # print(key)
      mapper[key] = data[[5]][i]
      
    }
    
  }
  return (mapper)

}

# data_path <- 'input/mapping_data.tsv'
# 
# mapper <- map_slseq_nissl(data_path)
# 
# mapper["Puck_210817_02_Seurat.RDS"]
