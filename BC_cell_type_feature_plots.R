#making violin plots of cell type marker genes

library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)

BC <- readRDS(file='BC_comb.rds')

feat <- c('COL1A1','CDH11','RUNX2','TOP2A','PCNA','MKI67','ACAN','COL2A1',
              'SOX9','CTSK','MMP9','IL7R','CD3D','NKG7','CD74','CD14','FCGR3A',
              'COL1A1','LUM','DCN','CDH11','ACTA2','RGS5','THY1','CXCL12',
              'SFRP2','CDH11','MME','MYLPF','MYL1','PECAM1','VWF','FBLN1','PTH1R',
          'IFITM5','ACP5','CD3D','CD8A','CD4','S100A8','FCN1','PLVAP','DMP1','IBSP',
          'SPI1','RSAD2','IRF8','GZMB','JCHAIN','TPSAB1',
          'CPA3','TPSB2','MS4A2','IFIT3','IFIT1','IFIT2')

DefaultAssay(BC) <- "RNA"
Idents(BC) <- BC$integrated_snn_res.0.45
#p <- DimPlot(BC, reduction = 'umap', label = TRUE)
#ggsave("umap_BC-combined_res045.pdf", plot = p)
BC.markers <- FindAllMarkers(BC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(BC.markers,file="BC_combined_res045_clustermarkers_RNA.csv")




p <- DimPlot(BC, reduction = 'umap', label = TRUE)
ggsave("umap_BC-combined_res045_RNA.pdf", plot = p)

for (i in 1:length(feat)){
  p <- FeaturePlot(BC, features = feat[i])
  ggsave(paste0(feat[i],"_RNA_umap.pdf"), plot = p)
  p <- RidgePlot(BC, features = feat[i])
  ggsave(paste0(feat[i],"_RNA_RidgePlot.pdf"), plot = p)
  p <- VlnPlot(BC, features = feat[i],pt.size=F)
  ggsave(paste0(feat[i],"_RNA_VlnPlot.pdf"), plot = p)
}

#saveRDS(BC, file = 'BC_comb.rds')