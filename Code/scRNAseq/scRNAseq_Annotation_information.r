library(tidyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(scCustomize)
merged_seurat <- readRDS('Cramp1.combined_CSS_scaled_annotated.rds')
Idents(merged_seurat) <- 'seurat_clusters'
new.cluster.ids <- c("Erythroblast", # 0
                     "Erythroid progenitor", # 1
                     "Myeloid progenitor", # 2
                     "Erythroid progenitor", # 3
                     "Erythroblast" , # 4
                     "Hepatocyte", # 5
                     "Hepatic stellate cell", # 6
                     "Liver endothelial cell", # 7
                     "Neutrophil", # 8
                     "Erythroid progenitor" , # ？？ 9 #经过xcr确定
                     "Hepatocyte", # 10
                     "Early erythrocyte", # 11
                     "Early erythrocyte", # 12
                     "TE_rich", # 13
                     "Kupffer cell", # 14
                     "Erythroblast", # 15
                     "Mast cell", # 16
                     "Neutrophil", # 17
                     "Megakaryocyte") # 18
names(new.cluster.ids) <- levels(merged_seurat)
merged_seurat <- RenameIdents(merged_seurat, new.cluster.ids)
merged_seurat@meta.data$celltype <- merged_seurat@active.ident
merged_seurat$celltype.group <- paste(merged_seurat$celltype,
                                      merged_seurat$orig.ident,sep = "_")
DimPlot(merged_seurat, reduction = "umap", label = TRUE, pt.size = 0.5)
