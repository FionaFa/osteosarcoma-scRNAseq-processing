#in the current clustering, pericytes and MSCs form one cluster. To differentiate the 
#two sub-populations, the resoution would have to be set to over 1.0, which would
#result in extreme over-clustering. Therefore, I will sub-cluster the pericyte/MSC 
#cluster and then transfer the annotation to the overall object

library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)

BC <- readRDS(file='BC_comb.rds')

names(BC@graphs)

Idents(BC) <- BC$integrated_snn_res.0.45
BC <- FindSubCluster(BC, "9",'integrated_snn',subcluster.name = 'peri.msc', resolution =0.45)

#p <- DimPlot(BC, reduction = 'umap', group.by = 'peri.msc',label = TRUE)
#ggsave("umap_peri-msc.pdf", plot = p)

saveRDS(BC, file = "BC_combined_2.rds")

DefaultAssay(BC) <- "RNA"
Idents(BC) <- BC$integrated_snn_res.0.45

BC_peri <- subset(x=BC, idents = c("9"))

BC_peri <- SetIdent(BC_peri, value = BC_peri@meta.data$peri.msc)


p <- DimPlot(BC_peri, reduction = 'umap', label = TRUE)
ggsave("umap_peri-msc_test.pdf", plot = p)


#BC.markers <- FindAllMarkers(BC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#BC.markers %>%
#  group_by(cluster) %>%
#  slice_max(n = 2, order_by = avg_log2FC)
#write.table(BC.markers,file="BC_clustermarkers.csv")

feat <- c('RGS5','ACTA2','MME','THY1','CXCL12','SFRP2')

for (i in 1:length(feat)){
  p <- FeaturePlot(BC_peri, features = feat[i])
  ggsave(paste0(feat[i],"_peri_umap.pdf"), plot = p)
  p <- RidgePlot(BC_peri, features = feat[i])
  ggsave(paste0(feat[i],"_peri_RidgePlot.pdf"), plot = p)
  p <- VlnPlot(BC_peri, features = feat[i],pt.size=F)
  ggsave(paste0(feat[i],"_peri_VlnPlot.pdf"), plot = p)
}
