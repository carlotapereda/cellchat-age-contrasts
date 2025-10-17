#Updated on 25 August 2025 by Carlota Pereda Serras
#This now includes CellChatDBv2 - checked with the new version (not reflected in code, just redownloaded package and explored the DB)
## LIBRARIES
library(CellChat)
library(ggplot2)
library(patchwork)
library(igraph)
library(Seurat)
options(stringsAsFactors = FALSE)


## CC OBJECT FUNCTION
createCCObjv2 <- function(sex, age, genotype, cellChatFileName) {
  start_time <- Sys.time()
  
  # Subset the data according to the metadata
  cell.use <- rownames(meta)[meta$sex == sex & meta$age == age & meta$genotype == genotype]
  
  # Prepare input data for CelChat analysis
  data.input.filt <- data.input[, cell.use]
  meta.filt <- meta[cell.use, ]
  
  # Create a CellChat (CC) object
  cellchat <- createCellChat(object = data.input.filt, meta = meta.filt, group.by = "cell_type_identity")
  
  # Add cell info into meta slot of the object (Optional)
  cellchat <- addMeta(cellchat, meta = meta.filt)
  cellchat <- setIdent(cellchat, ident.use = "cell_type_identity") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  # Set database
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing expression data for CCC analysis
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 6) # do parallel
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probability and infer network
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Extract the CC network as a data frame
  df.net <- subsetCommunication(cellchat)
  
  # Infer the CCC at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # Save the CellChat object
  saveRDS(cellchat, file = cellChatFileName)
  
  
  end_time <- Sys.time()
  end_time - start_time
  
  return(end_time - start_time)
  
  # return(cellchat)
}





# Run on all data

# ~~~~~~~~~~~~Define the possible values for each condition
sex_values <- c("Male", "Female")
age_values <- c("06Mo", "12Mo", "18Mo")
genotype_values <- c("E33", "E44")


# ~~~~~~~~~~~~Load data
load("/Users/carlotapereda/Library/CloudStorage/Box-Box/~~~~MyFolder/~~MyLab/R_YadongMouseAPOE4/data/Emouse.RData")
rm(Emouse.markers)
data.input <- GetAssayData(Emouse, assay = "RNA", slot = "data") 
# A dataframe with rownames containing cell mata data
meta <- Emouse@meta.data
# Remove the cell type called "remove" and Choroid Plexus
meta <- subset(meta, cell_type_identity != "remove")
meta <- subset(meta, cell_type_identity != "Choroid.Plexus")
rm(Emouse) 

# ~~~~~~~~~~~~ multisession
options(future.globals.maxSize = 4 * 1024^3)  # 4 GB
future::plan("multisession", workers = 6)

# ~~~~~~~~~~~~Iterate over each combination of conditions
for (sex in sex_values) {
  for (age in age_values) {
    for (genotype in genotype_values) {
      
      cellChatFileName <- paste0("CC_DBv2_", sex, "_", age, "_", genotype, ".rds")
      start_time <- Sys.time()
      print(paste0("~~~~starting ", cellChatFileName))
      createCCObjv2(sex, age, genotype, cellChatFileName)
      print(paste0("~~~~finished creating CC object for ", cellChatFileName))
      end_time <- Sys.time()
      end_time - start_time
    }
  }
}



