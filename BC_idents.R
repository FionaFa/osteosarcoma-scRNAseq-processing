library(Seurat)
library(ggplot2)
library(magrittr)

BC <- readRDS(file='BC_combined.rds')

#issue is that all patient have the same annotation --> we use the 
#original seurat object and cell ids to annotate the patient ids 
#in the combined object
#load all individual patient files
BC10 <- readRDS(file="BC10.rds")
BC11 <- readRDS(file="BC11.rds")
BC16 <- readRDS(file="BC16.rds")
BC17 <- readRDS(file="BC17.rds")
BC20 <- readRDS(file="BC2.rds")
BC21 <- readRDS(file="BC20.rds")
BC22 <- readRDS(file="BC21.rds")
BC2 <- readRDS(file="BC22.rds")
BC3 <- readRDS(file="BC3.rds")
BC5 <- readRDS(file="BC5.rds")
BC6 <- readRDS(file="BC6.rds")

#get list of cell ids in combined file
BC.ids <- colnames(BC)
BC.ids <- strsplit(BC.ids, '_')

#match patient ids with cell ids in combined file
pat.ids <- c()
for (i in 1:length(BC.ids)){
  print(i)
  if (BC.ids[[i]][1] %in% colnames(BC10)){
    pat.ids[i] <- 'BC10'
  } else if (BC.ids[[i]][1] %in% colnames(BC11)){
    pat.ids[i] <- 'BC11'
  } else if (BC.ids[[i]][1] %in% colnames(BC16)){
    pat.ids[i] <- 'BC16'
  } else if (BC.ids[[i]][1] %in% colnames(BC17)){
    pat.ids[i] <- 'BC17'
  } else if (BC.ids[[i]][1] %in% colnames(BC20)){
    pat.ids[i] <- 'BC20'
  } else if (BC.ids[[i]][1] %in% colnames(BC21)){
    pat.ids[i] <- 'BC21'
  } else if (BC.ids[[i]][1] %in% colnames(BC22)){
    pat.ids[i] <- 'BC22'
  } else if (BC.ids[[i]][1] %in% colnames(BC2)){
    pat.ids[i] <- 'BC2'
  } else if (BC.ids[[i]][1] %in% colnames(BC3)){
    pat.ids[i] <- 'BC3'
  } else if (BC.ids[[i]][1] %in% colnames(BC5)){
    pat.ids[i] <- 'BC5'
  } else if (BC.ids[[i]][1] %in% colnames(BC6)){
    pat.ids[i] <- 'BC6'
  }
}

write.csv(pat.ids, file = 'patids.csv')

#create metadata column with patient ids
BC <- AddMetaData(BC, metadata = pat.ids)
p <- DimPlot(BC, group.by = "orig.ident", reduction = "umap")
ggsave("umap_BC-combined_origident_2.pdf", plot = p)
saveRDS(BC, file="BC_comb_patids.rds")

