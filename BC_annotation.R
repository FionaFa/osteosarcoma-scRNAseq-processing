#annotation of cell types

library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)

BC <- readRDS(file='BC_combined_2.rds')

BC <- SetIdent(BC, value = BC@meta.data$peri.msc)

levels(BC)

#order of clusters: 2,6,3,7,10,1,12,13,8,4,0,9_0,9_1,9_2,9_3,11,14,9_5,17,15,16,9_4,5,9_6
celltypes <- c('chondroblastic cells','T cells','osteogenic cells','proliferating cells',
               'chondroblastic cells','fibroblasts','endothelial cells','macrophages',
               'proliferating cells','osteoclasts','myeloid cells','pericytes',
               'pericytes','pericytes','pericytes','peripheral monocyte cells',
               'inflamed osteoblastic cells','mesenchymal stem cells','mast cells',
               'dendritic cells','osteoblastic cells','mesenchymal stem cells',
               'osteocytes','pericytes')

names(celltypes) <- levels(BC)
BC <- RenameIdents(BC, celltypes)

p <- DimPlot(BC, reduction = 'umap',label = TRUE,label.size = 2) + NoLegend()
ggsave("umap_celltypes.pdf", plot = p)

BC$celltype <- Idents(BC)

saveRDS(BC, file = "BC_combined_3.rds")