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

# Figure S1a - Multiome dataset UMAP
# Color conversion functions
hsv2rgb <- function(x){  
  h <- x[1,1]  
  s <- x[2,1]  
  v <- x[3,1]    
  C <- s*v   
  hdash <- h*6  
  X <- C * (1 - abs(hdash %% 2 -1))
  if (0 <= hdash & hdash <=1) RGB1 <- c(C, X, 0)  
  if (1 <= hdash & hdash <=2) RGB1 <- c(X, C, 0)  
  if (2 <= hdash & hdash <=3) RGB1 <- c(0, C, X)  
  if (3 <= hdash & hdash <=4) RGB1 <- c(0, X, C)  
  if (4 <= hdash & hdash <=5) RGB1 <- c(X, 0, C)  
  if (5 <= hdash & hdash <=6) RGB1 <- c(C, 0, X)    
  RGB1 + (v-C)
}

pastellize <- function(x, p){
  if (is.character(x)) x <- col2rgb(x)/255
  if (is.numeric(x)) x <- matrix(x, nr=3)
  col <- rgb2hsv(x, maxColorValue=1)
  col[2,1] <- col[2,1]*p
  col <- hsv2rgb(col)
  rgb(col[1], col[2], col[3])
}

# Color definitions
purples <- pal_material("deep-purple", alpha = 0.5)(10)
blues <- pal_material("blue", alpha = 0.5)(10)
indigos <- pal_material("indigo", alpha = 0.5)(10)
reds <- pal_material("red", alpha = 0.5)(10)
greens <- pal_material("green", alpha = 0.5)(10)
oranges <- pal_material("orange", alpha = 0.5)(10)
blue_greys <- pal_material("blue-grey", alpha = 0.5)(10)
light_blues <- pal_material("light-blue", alpha = 0.5)(10)
greys <- pal_material("grey", alpha = 0.5)(10)
browns <- pal_material("brown", alpha = 0.5)(10)
pinks <- pal_material("pink", alpha = 0.5)(10)

colours_multiome_lvl1 <- c('TAL'=pastellize(indigos[5], 0.8),
                           'T Cell'=pastellize(browns[7], 0.5),
                           'Myeloid Cell'=pastellize(oranges[8], 0.5),
                           'Endothelia'=pastellize(reds[5], 0.8),
                           'PT'=pastellize(purples[5], 0.8),
                           'IC-B'=pastellize(greens[3], 0.8),
                           'IC-A'=pastellize(greens[7], 0.8),
                           'CNT'=pastellize(blues[5], 0.8),
                           'B Cell'=pastellize(pinks[5], 0.3),
                           'DCT'=pastellize(blues[9], 0.8),
                           'PC'=pastellize(blues[2], 0.8),
                           'Interstitium'=pastellize(light_blues[3], 0.8),
                           'PEC'=pastellize(greys[6], 0.3),
                           'ATL'=pastellize(blue_greys[2], 0.8),
                           'Podocyte'=pastellize(greys[8], 0.3),
                           'DTL'=pastellize(blue_greys[4], 0.8))

p <- DimPlot(multiome, label=F, pt.size=0.1, cols=colours_multiome_lvl1, group.by = 'Annotation.Lvl1', reduction='umap_wnn', order=F) + NoLegend() + NoAxes() + ggtitle('')
LabelClusters(p, id = "Annotation.Lvl1", size=5, fontface = "bold", color = "black", box=F, repel=F)

ggsave(filename = file.path(path, 'multiome_umap.png'),
       scale = 0.5, width = 36, height = 30, units='cm')


# Figure 1b - Proximal tubule markers
FeaturePlot(multiome_pt, "CUBN", reduction='DM', order=T, cols=c('grey90', 'navy')) + ggtitle('') + ggtitle('') + NoAxes()
ggsave(filename = file.path(path, 'umap_cubn.png'),
       scale = 0.5, width = 25, height = 20, units='cm')

FeaturePlot(multiome_pt, "MME", reduction='DM', order=T, cols=c('grey90', 'navy')) + ggtitle('') + ggtitle('') + NoAxes()
ggsave(filename = file.path(path, 'umap_mme.png'),
       scale = 0.5, width = 25, height = 20, units='cm')

FeaturePlot(multiome_pt, "HAVCR1", reduction='DM', order=T, cols=c('grey90', 'navy')) + ggtitle('') + ggtitle('') + NoAxes()
ggsave(filename = file.path(path, 'umap_havcr1.png'),
       scale = 0.5, width = 25, height = 20, units='cm')

FeaturePlot(multiome_pt, "VCAM1", reduction='DM', order=T, cols=c('grey90', 'navy')) + ggtitle('') + ggtitle('') + NoAxes()
ggsave(filename = file.path(path, 'umap_vcam1.png'),
       scale = 0.5, width = 25, height = 20, units='cm')

FeaturePlot(multiome_pt, "CCL2", reduction='DM', order=T, cols=c('grey90', 'navy')) + ggtitle('') + ggtitle('') + NoAxes()
ggsave(filename = file.path(path, 'umap_ccl2.png'),
       scale = 0.5, width = 25, height = 20, units='cm')

FeaturePlot(multiome_pt, "CDKN1A", reduction='DM', order=T, cols=c('grey90', 'navy')) + ggtitle('') + ggtitle('') + NoAxes()
ggsave(filename = file.path(path, 'umap_cdkn1a.png'),
       scale = 0.5, width = 25, height = 20, units='cm')


# Figure 1c - RPTEC irradiation volcano plot
DEGs <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
DEGs <- as.data.frame(DEGs)
DEGs$gene_id <- rownames(DEGs)

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                 values = rownames(DEGs), mart = ensembl) 
DEGs <- merge(DEGs, genemap, by.x="gene_id", by.y="ensembl_gene_id")
DEGs$hgnc_symbol <- make.unique(DEGs$hgnc_symbol)
DEGs$log2FoldChange <- DEGs$log2FoldChange * (-1)

EnhancedVolcano(DEGs,
                lab = DEGs$hgnc_symbol,
                selectLab = c(''),
                x = 'log2FoldChange',
                y = 'pvalue',
                legendPosition = 'none',
                pCutoff = 0.05,
                FCcutoff = 1,
                subtitle='',
                caption='', drawConnectors=T, lengthConnectors = unit(2, "npc"),
                widthConnectors = 1.0,
                colConnectors = 'black',
                arrowheads = FALSE,
                col=c('grey60', 'grey60', 'grey60', 'red4'), colAlpha=1) + ggtitle('') + 
  theme(axis.text.x = element_text(size = 10, color='black'),
        axis.text.y = element_text(size = 10, color='black'))

ggsave(filename = file.path(path, 'rptec_volcano.png'),
       scale = 0.5, width = 20, height = 20, units='cm')
