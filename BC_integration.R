library(utils)
library(Seurat)
library(magrittr)
library(dplyr)
library(DoubletFinder)
library(parallel)
library(utils)
library(infercnv)
library(R.utils)
library(readr)
library(ggplot2)

#load datasets
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

#make a list object with all the data
BC <- list(BC10, BC11, BC16, BC17, BC20, BC21, BC22, BC2, BC3, BC5, BC6)
names(BC) <- c('BC10', 'BC11', 'BC16', 'BC17', 'BC20', 'BC21', 'BC22', 'BC2', 'BC3', 'BC5', 'BC6')
#doubletremoval
for (i in 1:length(BC)) {
  BCt <- BC[[i]]
  BCt <- RunPCA(BCt, features = VariableFeatures(object = BCt))
  #find significant pcs
  stdv <- BCt[["pca"]]@stdev
  sum.stdv <- sum(BCt[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  BCt <- FindNeighbors(BCt, dims = 1:min.pc)
  BCt <- FindClusters(BCt, resolution = 0.45)
  BCt <- RunUMAP(BCt, dims = 1:min.pc)
  p <- DimPlot(BCt, reduction = "umap")
  ggsave(paste0(names(BC)[i],"_umap_before_doubletremoval.pdf"))
  #doublet identification and removal
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(BCt, PCs = 1:min.pc)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- BCt@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(BCt@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  BCt <- doubletFinder_v3(seu = BCt, 
                           PCs = 1:min.pc, 
                           pK = optimal.pk,
                           nExp = nExp.poi.adj)
  metadata <- BCt@meta.data
  colnames(metadata)[8] <- "doublet_finder"
  BCt@meta.data <- metadata 
  
  # subset and save
  BCt.singlets <- subset(BCt, doublet_finder == "Singlet")
  BCt <- BCt.singlets
  remove(BCt.singlets)
  
  #umap after doublet removal
  p2 <- DimPlot(BCt, reduction = "umap")
  ggsave(paste0(names(BC)[i], "umap_after_doubletremoval.pdf"))
  BC[i] <- BCt
}

#Integration
features <- SelectIntegrationFeatures(object.list = BC)

BC.anchors <- FindIntegrationAnchors(object.list = BC, anchor.features = features)

BC.combined <- IntegrateData(anchorset = BC.anchors)

saveRDS(BC.combined, file = "BC_combined.rds")

#processing
DefaultAssay(BC.combined) <- "integrated"

#clustering
BC.combined <- ScaleData(BC.combined, verbose = FALSE)
BC.combined <- RunPCA(BC.combined, npcs = 30, verbose = FALSE)

##find number of significant pcs
stdv <- BC.combined[["pca"]]@stdev
sum.stdv <- sum(BC.combined[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)


BC.combined <- RunUMAP(BC.combined, reduction = 'pca', dims = 1:min.pc)
BC.combined <- FindNeighbors(BC.combined, reduction = 'pca', dims = 1:min.pc)
BC.combined <- FindClusters(BC.combined, resolution = 0.45)

p <- DimPlot(BC.combined, reduction = 'umap', label = TRUE)

ggsave("umap_BC-combined.png", plot = p)

BC.markers <- FindAllMarkers(BC.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.table(BC.markers,file="BC_combined_clustermarkers.csv")

saveRDS(BC.combined, file="BC_combined.rds")

