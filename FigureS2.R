library(Seurat)
library(Signac)
library(SeuratObject)
library(SeuratDisk)
library(SeuratWrappers)
library(scCustomize)
library(UCell)
library(harmony)
library(CellChat)
library(destiny)
library(DESeq2)
library(clusterProfiler)
library(EnhancedVolcano)
library(ggplot2)
library(ggsci)
library(destiny)
library(scales)
library(viridis)
library(nichenetr)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)

multiome <- readRDS('/Users/maximilianreck/Drive/Datasets/seurat_object_multiome.rds')
#-------------------------------------------------------------------------------

# Figure S2a - Heatmap of macrophage markers
# Average expression across clusters
Idents(multiome_mac) <- multiome_mac$Cluster
average_expr <- AverageExpression(multiome_mac)
average_expr <- as.data.frame(average_expr$SCT)

average_expr <- average_expr[rownames(average_expr)%in%c('IL6' 'CD36', 'LYVE1', 'RNASE1', 'NFKBIZ', 'CCL8', 'F13A1', 'CCL2', 'CD80',
                                                         'MAF', 'CD209', 'MRC1', 'CXCL10', 'CCL4', 'CD274', 'CXCL9', 'IRF1', 'NFKB1', 'CCL21', 
                                                         'ICAM1', 'CCL3', 'CD44', 'CD58', 'VCAM1', 'IL1B', 'CD83', 'IRF8', 'APOE', 'SELENOP', 
                                                         'CD163', 'CD38', 'MMP9', 'CXCL12', 'PLAU', 'C1QA', 'MERTK', 'CD206'),]

# Plot
pheatmap(t(average_expr), cluster_rows=T, cluster_cols=T, show_rownames=T, show_colnames=T,
         color=cols, 
         clustering_method="ward.D2",
         scale="column", 
         fontsize_row=12, 
         labels_col = parse(text = paste0("italic('", rownames(average_expr), "')")))
         gaps_col = 5, gaps_row = 9)

