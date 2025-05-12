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

# Figure 1a - Proximal tubule diffusion map
# Extract PT cells and format
multiome_pt <- subset(multiome, subset=Annotation.Lvl1 %in% c('PT'))
dr <- multiome_pt@reductions[["umap_wnn"]]
multiome_pt <- CreateSeuratObject(counts = multiome_pt@assays$RNA@counts, meta.data = multiome_pt@meta.data)
multiome_pt@reductions[["umap"]] <- dr
multiome_pt <- SCTransform(multiome_pt, vst.flavor = "v2")
multiome_pt <- RunPCA(multiome_pt, npcs = 100)

multiome_pt$celltype <- multiome_pt$Annotation.Lvl2
multiome_pt$celltype[multiome_pt$celltype=='PT S1'] <- 'S1'
multiome_pt$celltype[multiome_pt$celltype=='PT S2'] <- 'S2'
multiome_pt$celltype[multiome_pt$celltype=='PT S3'] <- 'S3'
multiome_pt$celltype[multiome_pt$celltype=='PT Injured'] <- 'Altered'
multiome_pt$celltype[multiome_pt$celltype=='PT Inflammatory'] <- 'Altered'

indigos <- pal_material("indigo", alpha = 1)(10)
oranges <- pal_material("orange", alpha = 1)(10)


# Calculate diffusion map
pca <- as.data.frame(t(multiome_pt@reductions[["pca"]]@cell.embeddings))
colnames(pca) <- multiome_pt$Annotation.Lvl2

dm <- DiffusionMap(t(pca[1:50,]))
dpt <- DPT(dm, tips=c(2746))

# Add embedding to Seurat object
dr <- multiome_pt@reductions[["umap"]]
dr@cell.embeddings[,1] <- as.numeric(dm$DC1)*100
dr@cell.embeddings[,2] <- as.numeric(dm$DC2)*100
dr@key <- 'DC_'
colnames(dr@cell.embeddings) <- c('DC_1', 'DC_2')
multiome_pt@reductions[["DM"]] <- dr

# Add pseudotime values
multiome_pt$pseudotime <- (rank(dpt$dpt) / length(dpt$dpt))
multiome_pt$pseudotime <- (-multiome_pt$pseudotime)+1


# Plot
DimPlot(multiome_pt, reduction='DM', group.by = 'celltype',
        cols = c('S1'=indigos[5], 'S2'=indigos[3], 'S3'=indigos[8], 'Altered'=oranges[4])) + 
  ggtitle('') + NoAxes() + theme_gray()

ggsave(filename = file.path(path, 'umap_pt.pdf'),
       scale = 0.5, width = 25, height = 20, units='cm')


# Figure 1b - Diffusion map by condition
# Format meta data
multiome_pt$condition_plot <- multiome_pt$Condition
multiome_pt$condition_plot[multiome_pt$condition_plot=='UUO'] <- 'Injured'
multiome_pt$condition_plot[multiome_pt$condition_plot=='Control'] <- 'Normal'
multiome_pt$condition_plot <- factor(multiome_pt$condition_plot, levels=c('Normal', 'Injured'))

# Plot
DimPlot(multiome_pt, reduction='DM', group.by = 'condition_plot',
        cols = c('grey70', 'navy'), order=T) + 
  ggtitle('') + NoAxes() + theme_void()

ggsave(filename = file.path(path, 'umap_condition.pdf'),
       scale = 0.5, width = 25, height = 20, units='cm')


# Figure 1c - Diffusion map by pseudotime
# Plot
FeaturePlot(multiome_pt, features = 'pseudotime', reduction='DM', order=F) + scale_color_viridis() + ggtitle('') + NoAxes()

ggsave(filename = file.path(path, 'umap_pseudotime.pdf'),
       scale = 0.5, width = 25, height = 20, units='cm')


# Figure 1d - Senescence gene set scores
# Load gene sets
sen_mayo <- c('ACVR1B', 'ANG', 'ANGPT1', 'ANGPTL4', 'AREG', 'AXL', 'BEX3', 'BMP2', 'BMP6', 'C3', 'CCL1', 'CCL13', 'CCL16', 'CCL2', 'CCL20', 'CCL24', 'CCL26', 'CCL3', 'CCL3L1', 'CCL4', 'CCL5', 'CCL7', 'CCL8', 'CD55', 'CD9', 'CSF1', 'CSF2', 'CSF2RB', 'CST4', 'CTNNB1', 'CTSB', 'CXCL1', 'CXCL10', 'CXCL12', 'CXCL16', 'CXCL2', 'CXCL3', 'CXCL8', 'CXCR2', 'DKK1', 'EDN1', 'EGF', 'EGFR', 'EREG', 'ESM1', 'ETS2', 'FAS', 'FGF1', 'FGF2', 'FGF7', 'GDF15', 'GEM', 'GMFG', 'HGF', 'HMGB1', 'ICAM1', 'ICAM3', 'IGF1', 'IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5', 'IGFBP6', 'IGFBP7', 'IL10', 'IL13', 'IL15', 'IL18', 'IL1A', 'IL1B', 'IL2', 'IL32', 'IL6', 'IL6ST', 'IL7', 'INHA', 'IQGAP2', 'ITGA2', 'ITPKA', 'JUN', 'KITLG', 'LCP1', 'MIF', 'MMP1', 'MMP10', 'MMP12', 'MMP13', 'MMP14', 'MMP2', 'MMP3', 'MMP9', 'NAP1L4', 'NRG1', 'PAPPA', 'PECAM1', 'PGF', 'PIGF', 'PLAT', 'PLAU', 'PLAUR', 'PTBP1', 'PTGER2', 'PTGES', 'RPS6KA5', 'SCAMP4', 'SELPLG', 'SEMA3F', 'SERPINB4', 'SERPINE1', 'SERPINE2', 'SPP1', 'SPX', 'TIMP2', 'TNF', 'TNFRSF10C', 'TNFRSF11B', 'TNFRSF1A', 'TNFRSF1B', 'TUBGCP2', 'VEGFA', 'VEGFC', 'VGF', 'WNT16', 'WNT2')
tabula_muris <- c("0910001L09Rik", "2010107E04Rik", "2310044H10Rik", "2900097C17Rik", "9530068E07Rik", "Aco2", "Actb", "Actg1", "Actr3", "Adam10", "Add1", "Adrm1", "Aes", "Aff4", "Amfr", "Anapc11", "Anxa2", "Ap2s1", "Aplp2", "Apoe", "App", "Arcn1", "Arhgdia", "Arpc1b", "Arpc3", "Asah1", "Atp1a1", "Atp1b3", "Atp5b", "Atp5e", "B2m", "Bri3", "Bsg", "Btg1", "Btg2", "Bzw1", "Cab39", "Calm2", "Calr", "Calu", "Canx", "Capn2", "Caprin1", "Capza1", "Capza2", "Cat", "Ccni", "Ccpg1", "Cd164", "Cd63", "Cd74", "Cd9", "Cdc42", "Cfl1", "Chmp2a", "Cirbp", "Clic1", "Clic4", "Clk1", "Clptm1l", "Cltc", "Cmtm6", "Copa", "Copb1", "Cotl1", "Cox6a1", "Cox6b1", "Cox6c", "Cox8a", "Csde1", "Cst3", "Ctnnb1", "Ctsd", "Cyba", "Cyfip1", "Dcn", "Ddx3x", "Ddx5", "Ddx6", "Dnaja2", "Dnm2", "Dpysl2", "Drap1", "Dusp1", "Eef1a1", "Eef2", "Egr1", "Eif3a", "Eif3b", "Eif3f", "Eif4a2", "Eif4g2", "Eif5a", "Erbb2ip", "Errfi1", "Fau", "Fos", "Fth1", "Ftl1", "Gadd45gip1", "Galnt1", "Gapdh", "Gas6", "Gdi2", "Ghitm", "Glg1", "Glud1", "Glul", "Gm1821", "Gnai3", "Gnb1", "Gnb2l1", "Golph3", "Gpbp1", "Gpx3", "Gpx4", "Grn", "Gsk3b", "Gsn", "Gstm1", "H2-Aa", "H2-Ab1", "H2-D1", "H2-K1", "Hbp1", "Hif1a", "Hist1h4d", "Hnrnpa2b1", "Hnrnpf", "Hnrnph1", "Hnrnpk", "Hnrnpu", "Hsd17b4", "Hspa5", "Hspa8", "Hspa9", "Ier2", "Ifitm2", "Iqgap1", "Itfg1", "Itgb1", "Ivns1abp", "Jmjd1c", "Jun", "Junb", "Jund", "Kdm6b", "Klf13", "Lamp1", "Lamp2", "Laptm4a", "Lars2", "Lrrc58", "Lrrc8a", "Lyz2", "M6pr", "Macf1", "Malat1", "Map1lc3a", "Map3k1", "Mapk1", "Mat2a", "Matr3", "Mbnl1", "Mcl1", "Msn", "Mtpn", "Myadm", "Ndufa1", "Nedd4", "Nfic", "Nono", "Npm1", "Nptn", "Nucb1", "Nufip2", "Oaz1", "Ogdh", "Osbpl9", "P4hb", "Pan3", "Pbrm1", "Pdcd4", "Pdcd6ip", "Pdia3", "Pdia4", "Pebp1", "Pfn1", "Ppia", "Ppp2ca", "Ppp2r1a", "Prkar1a", "Pros1", "Prr13", "Prrc2c", "Psap", "Psmd2", "Psmd4", "Ptma", "Ptms", "Ptpra", "Pttg1ip", "Rab18", "Rab6a", "Rabac1", "Rac1", "Rap1b", "Reep5", "Rhoa", "Rhoc", "Rlim", "Rnaset2b", "Rod1", "Rpl13", "Rpl13a", "Rpl14", "Rpl17", "Rpl23", "Rpl23a", "Rpl3", "Rpl31", "Rpl32", "Rpl36", "Rpl37", "Rpl37a", "Rpl39", "Rpl4", "Rpl41", "Rpl6", "Rpl9", "Rplp0", "Rplp1", "Rpn1", "Rpn2", "Rps10", "Rps12", "Rps14", "Rps16", "Rps17", "Rps18", "Rps19", "Rps20", "Rps21", "Rps23", "Rps25", "Rps26", "Rps27", "Rps28", "Rps29", "Rps3", "Rps3a", "Rps5", "Rps7", "Rpsa", "Rrp1", "Sdc4", "Sdcbp", "Sdha", "Sdhc", "Sec61a1", "Sepp1", "Sept2", "Sepw1", "Serbp1", "Serinc1", "Serinc3", "Serp1", "Sf3b2", "Sf3b4", "Sh3bgrl", "Shfm1", "Shisa5", "Skp1a", "Slc39a1", "Slc44a2", "Snrpc", "Sod1", "Sparc", "Spnb2", "Spr", "Srsf5", "Ssbp4", "Ssr1", "Ssr3", "Stt3a", "Sun2", "Surf4", "Swi5", "Sypl", "Tagln2", "Taok1", "Tax1bp1", "Tex261", "Tm9sf2", "Tm9sf3", "Tmbim6", "Tmed9", "Tmem123", "Tmem160", "Tmem30a", "Tmsb10", "Tmsb4x", "Tnks2", "Tomm6", "Tomm7", "Top2b", "Tpm4", "Tpp1", "Tpt1", "Tram1", "Trf", "Tsc22d3", "Txnip", "Uba1", "Uba2", "Uba52", "Ubb", "Ubc", "Ube2m", "Ubxn1", "Usmg5", "Vcp", "Vps35", "Wbp2", "Wdr1", "Ybx1", "Ywhaz", "Zfp106", "Zfp207", "Zfp36")
tabula_muris <- convert_mouse_to_human_symbols(tabula_muris)
tabula_muris <- tabula_muris[!is.na(tabula_muris)]

# Add module scores
multiome_pt <- AddModuleScore_UCell(multiome_pt, features = list('sen_mayo'= sen_mayo),
                                  chunk.size = 8000, ncores = 10, name='')
multiome_pt <- AddModuleScore_UCell(multiome_pt, features = list('tabula_muris'= tabula_muris),
                                  chunk.size = 8000, ncores = 10, name='')

# Group by pseudotime
multiome_pt$class <- multiome_pt$pseudotime
multiome_pt$class[multiome_pt$class>0.85] <- 'PT senescence'
multiome_pt$class[multiome_pt$class>0.5 & multiome_pt$class<=0.85] <- 'PT acute injury'
multiome_pt$class[multiome_pt$class<=0.5] <- 'PT healthy'
multiome_pt$class <- factor(multiome_pt$class, levels=rev(c('PT healthy', 'PT acute injury', 'PT senescence')))

# Plot
DotPlot(multiome_pt, features=c('tabula_muris', 'sen_mayo'), group.by = 'class', 
        scale=T, col.min=-2, scale.min=100) + scale_color_viridis() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) + 
  RotatedAxis() + scale_x_discrete(labels= c('Tabula Muris Senis', 'SenMayo'))

ggsave(filename = file.path(path, 'signature_scores.pdf'),
       scale = 0.5, width = 19, height = 15, units='cm')


# Figure 1f - Heatmap irradiated RPTEC
# Load counts and run DESeq
counts <- read.csv(file.path(path, 'rptec_irr_counts.csv'))
rownames(counts) <- counts$X; counts$X <- NULL
counts <- counts[rowSums(counts)>0,]
rownames(counts) <- make.unique(gsub("\\..*","",rownames(counts)))
meta <- data.frame(names = c('B1-NR', 'B1-IR', 'B3-NR', 'B3-IR', 'B4-NR', 'B4-IR', 'B5-NR', 'B5-IR', 'B6-NR', 'B6-IR'),
                   Replicate = c('B1', 'B1', 'B3', 'B3', 'B4', 'B4', 'B5', 'B5', 'B6', 'B6'),
                   Condition = c('NR', 'IR', 'NR', 'IR', 'NR', 'IR', 'NR', 'IR', 'NR', 'IR'))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ Condition+Replicate)
dds <- dds[rowSums(counts(dds)) >= 100,]
dds <- DESeq(dds)
normalized_counts <- DESeq2::vst(dds, blind=FALSE)

# Translate gene IDs
rptec_matrix <- as.data.frame(assay(normalized_counts))
rptec_matrix$gene_id <- rownames(rptec_matrix)

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                 values = rownames(rptec_matrix), mart = ensembl) 
rptec_matrix <- merge(rptec_matrix, genemap, by.x="gene_id", by.y="ensembl_gene_id")
rownames(rptec_matrix) <- make.unique(rptec_matrix$hgnc_symbol)
rptec_matrix$hgnc_symbol <- NULL
rptec_matrix$gene_id <- NULL
colnames(rptec_matrix) <- sub('\\.', '-', colnames(rptec_matrix))

# Plot
genes <- unique(c('CDKN1A', 'MKI67', 'CCL2', 'CCL28', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 'CXCL8', 'CXCL16'))
rptec_matrix_ss <- rptec_matrix[(rownames(rptec_matrix) %in% genes),]

annotations <- as.data.frame(rep(c('Control', 'Irradiated'), 5))
colnames(annotations) <- c('Treatment')
row.names(annotations) <- colnames(rptec_matrix_ss)

annot_colors=list(Treatment=c(Irradiated="purple4", Control="grey70"))
rptec_matrix_ss <- rptec_matrix_ss[match(genes, rownames(rptec_matrix_ss)),]
rptec_matrix_ss <- rptec_matrix_ss[,match(c('B1-NR', 'B3-NR', 'B4-NR', 'B5-NR', 'B6-NR',
                                            'B1-IR', 'B3-IR', 'B4-IR', 'B5-IR', 'B6-IR'), colnames(rptec_matrix_ss))]

pheatmap(rptec_matrix_ss, 
         scale = 'row', 
         color = colorRampPalette(c(muted("navy", l=30, c=70), "white", muted("red", l=40, c=90)))(500),
         annotation_col = annotations, 
         annotation_colors = annot_colors, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         clustering_method = 'ward.D2', 
         gaps_row = 2, 
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         fontsize = 12,
         labels_row = parse(text = paste0("italic('", rownames(rptec_matrix_ss), "')")),
         labels_col = colnames(rptec_matrix_ss),
         gaps_col = c(5), 
         border_color = 'grey50')


# Figure 1g - Dotplot irradiated RPTEC markers in snRNA-seq
multiome_pt$class <- multiome_pt$pseudotime
multiome_pt$class[multiome_pt$class>0.85] <- 'PT senescence'
multiome_pt$class[multiome_pt$class>0.5 & multiome_pt$class<=0.85] <- 'PT acute injury'
multiome_pt$class[multiome_pt$class<=0.5] <- 'PT healthy'
multiome_pt$class <- factor(multiome_pt$class, levels=(c('PT healthy', 'PT acute injury', 'PT senescence')))

DotPlot(multiome_pt, 
        features = rev(c('CDKN1A', 'CCL2', 'CCL20', 'CCL28', 
                         'CXCL1', 'CXCL2', 'CXCL3', 'CXCL6', 
                         'CXCL8', 'CXCL12', 'CXCL16')), 
        group.by = 'class', 
        cols = c('grey70', 'navy')) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic")) +
  RotatedAxis() + 
  NoLegend() + 
  coord_flip()

ggsave(filename = file.path(path, 'senescence_markers_dotplot.pdf'),
       scale = 0.5, width = 10, height = 20, units='cm')


# Figure 1h - RPTEC genes vs pseudotime
# Find DEGs
DEGs <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
DEGs <- as.data.frame(DEGs)
DEGs$gene_id <- rownames(DEGs)

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                 values = rownames(DEGs), mart = ensembl) 
DEGs <- merge(DEGs, genemap, by.x="gene_id", by.y="ensembl_gene_id")

DEGs$hgnc_symbol <- make.unique(DEGs$hgnc_symbol)

DEGs$log2FoldChange <- DEGs$log2FoldChange * (-1)
DEGs <- DEGs[DEGs$hgnc_symbol %in% rownames(multiome_pt),]
DEGs <- DEGs[order(DEGs$log2FoldChange, decreasing = T),]
genes_up <- DEGs[DEGs$log2FoldChange>(1) & DEGs$padj < 0.05, 'hgnc_symbol']

# Add module score
multiome_pt <- AddModuleScore_UCell(multiome_pt,features = list('upregulated'=genes_up),
                                  chunk.size = 8000, ncores = 10, name='')
multiome_pt$upregulated_scaled <- scale(multiome_pt$upregulated)

matrix <- as.data.frame(t(GetAssayData(object = multiome_pt, assay = "SCT", slot = "data")))
matrix$pseudotime <- multiome_pt$pseudotime
matrix$rptec <- multiome_pt$upregulated_scaled

# Plot
ggplot(matrix) +
  geom_smooth(data=matrix, aes(x=pseudotime, y=rptec),
              method = "loess", se=F, span=0.7, color='cornflowerblue', size=1.5) +
  xlab('') + ylab('') +
  theme(axis.text.x = element_text(size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        legend.title = element_text(size=12, color="black"),
        legend.text = element_text(size=12, color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

ggsave(filename = file.path(path, 'rptec_score.pdf'),
       scale = 0.5, width = 15, height = 10, units='cm')


# Figure 1l - Myeloid cells UMAP
# Subset and calculate UMAP
multiome_myeloid <- subset(multiome, subset=Annotation.Lvl1 %in% c('Myeloid Cell'))
multiome_myeloid <- subset(multiome_myeloid, subset=Annotation.Lvl2 %in% c('pDC', 'Mast Cell'), invert=T)

multiome_myeloid <- SCTransform(multiome_myeloid, vst.flavor = "v2")
multiome_myeloid <- RunPCA(multiome_myeloid, npcs = 100)
multiome_myeloid <- RunHarmony(multiome_myeloid, assay.use="SCT", group.by.vars = "Library")
multiome_myeloid <- RunUMAP(multiome_myeloid, reduction = "harmony", dims = 1:50, min.dist = 0.2, n.neighbors=50)

# Remove small artifact cluster
select.cells <- CellSelector(plot = DimPlot(multiome_myeloid, label=T, reduction='umap'))
Idents(multiome_myeloid, cells = select.cells) <- "selected"
multiome_myeloid <- subset(multiome_myeloid, idents='selected', invert=T)
multiome_myeloid <- RunUMAP(multiome_myeloid, reduction = "harmony", dims = 1:50, min.dist = 0.2, n.neighbors=40)

# Cluster
multiome_myeloid <- FindNeighbors(multiome_myeloid, reduction = "harmony", dims = 1:50) %>% FindClusters(algorithm=3, resolution = 2)
Idents(multiome_myeloid) <- paste('c', Idents(multiome_myeloid), sep='')

multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c8` = "7")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c13` = "7")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c6` = "8")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c5` = "8")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c10` = "9")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c15` = "9")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c0` = "1")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c1` = "2")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c3` = "3")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c9` = "3")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c12` = "4")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c2` = "5")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c17` = "5")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c11` = "6")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c7` = "10")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c4` = "10")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c14` = "11")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c16` = "12")
multiome_myeloid$Cluster <- as.character(Idents(multiome_myeloid))

# Annotate cell types
Idents(multiome_myeloid) <- paste('c', Idents(multiome_myeloid), sep='')
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c8` = "CD16 Monocyte")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c13` = "CD16 Monocyte")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c6` = "CD14 Monocyte")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c5` = "Intermediate")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c10` = "Intermediate")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c15` = "Intermediate")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c11` = "Intermediate")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c2` = "Macrophage")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c17` = "Macrophage")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c0` = "Macrophage")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c1` = "Macrophage")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c3` = "Macrophage")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c9` = "Macrophage")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c12` = "Macrophage")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c7` = "cDC")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c4` = "cDC")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c14` = "cDC")
multiome_myeloid <- RenameIdents(object = multiome_myeloid, `c16` = "cDC")
multiome_myeloid$CellType <- as.character(Idents(multiome_myeloid))

# Color definitions
purples <- pal_material("deep-purple", alpha = 1)(10)
blues <- pal_material("blue", alpha = 1)(10)
indigos <- pal_material("indigo", alpha = 1)(10)
reds <- pal_material("red", alpha = 1)(10)
greens <- pal_material("green", alpha = 1)(10)
oranges <- pal_material("orange", alpha = 1)(10)
blue_greys <- pal_material("blue-grey", alpha = 1)(10)
light_blues <- pal_material("light-blue", alpha = 1)(10)
greys <- pal_material("grey", alpha = 1)(10)
browns <- pal_material("brown", alpha = 1)(10)
pinks <- pal_material("pink", alpha = 1)(10)

# Plot
multiome_myeloid$CellType <- factor(multiome_myeloid$CellType, levels=c('CD16 Monocyte', 'CD14 Monocyte',
                                                                    'Intermediate', 'cDC', 'Macrophage'))
DimPlot(multiome_myeloid, pt.size=0.8, order=F, group.by = 'CellType', reduction='umap',
        label=F, label.size=6, repel=T,
        cols=c('Macrophage'=purples[4],'CD14 Monocyte'=oranges[5], 'CD16 Monocyte'=oranges[8], 'Intermediate'=oranges[2], 
               'cDC'=browns[6])) + ggtitle('') + NoAxes() + theme_gray() #+ NoLegend()

ggsave(filename = file.path(path, 'umap_ct.pdf'),
       scale = 0.5, width = 30, height = 20, units='cm')


# Figure 1j - Dotplot of myeloid markers
multiome_myeloid$CellType <- factor(multiome_myeloid$CellType, levels=rev(c('CD16 Monocyte', 'CD14 Monocyte',
                                                                        'Intermediate', 'cDC', 'Macrophage')))
DotPlot(multiome_myeloid, 
        features = c('FCGR3A', 'VCAN', 'FCN1', 'SORL1', 'ITPR1', 'FLT3', 'IDO1', 'MRC1', 'CD163'),
        scale.min = 0, 
        group.by = 'CellType', 
        cols = c('grey90', 'navy')) +
  RotatedAxis() +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "white"))) +
  theme(axis.text.x = element_text(face = "italic"))

ggsave(filename = file.path(path, 'dotplot_myeloid.pdf'),
       scale = 0.5, width = 30, height = 20, units='cm')


DotPlot(subset_macrophage, features = c('CABLES1', 'CDKN1C', 'PDK4', 'ADAMTS15', 'MS4A6A', 'TGFB1', 'CEACAM4', 'THBS1', 'THBS1-IT1', 'CCL1', 'CCL3L3', 
                                        'IL12B', 'IL17A', 'CSF2', 'IL17F', 'CCL2', 'CCL4', 'CCL8', 'IL6', 'CLEC5A', 'TNFSF15', 'TNF', 'EDN1', 'HSD11B1', 'SERPINE1', 
                                        'SIMALR'),
        scale.min=0, group.by = 'Cluster', cols = c('grey90', 'navy')) +
  RotatedAxis() +
  guides(size=guide_legend(override.aes=list(shape=21, fill="white"))) 


# Figure 1k - Cellchat analysis PT -> Myeloid
# Format data
myeloid <- CreateSeuratObject(counts = multiome_myeloid@assays$RNA@counts, meta.data = multiome_myeloid@meta.data)
pt$CellType <- pt$class
pt <- CreateSeuratObject(counts = multiome_pt@assays$RNA@counts, meta.data = multiome_pt@meta.data)

data <- subset(multiome, cells=c(colnames(myeloid), colnames(pt)))
celltype <- data.frame(CellID=c(myeloid$CellID, pt$CellID), CellType=c(myeloid$CellType, pt$CellType))
celltype <- celltype[base::match(data$CellID, celltype$CellID),]
celltype$CellID==data$CellID
data$CellType <- celltype[base::match(data$CellID, celltype$CellID), 'CellType']

# Run cellchat
cellchat <- createCellChat(object = data@assays[["SCT"]]@data, 
                           meta = data@meta.data, 
                           group.by = 'CellType')

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 10)
options(future.globals.maxSize = 6000 * 1024^2)

cellchat <- identifyOverExpressedGenes(cellchat, min.cells=3)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)

computeAveExpr(cellchat, features = c("TNF","CXCL1"), type =  "truncatedMean", trim = 0.005)
cellchat <- computeCommunProb(cellchat, raw.use=F, population.size=F, type = "truncatedMean", trim = 0.005)


cellchat <- filterCommunication(cellchat, min.cells = 2)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Plot
netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(3, 2, 6, 4, 7, 3), thresh=0.01,
                 pairLR.use = data.frame(interaction_name=(c('CCL2_CCR2', 'CCL20_CCR6', 'CXCL12_CXCR4',
                                                             'ICAM1_ITGAX_ITGB2', 'VCAM1_INTEGRIN_ADB2', 
                                                             'TNF_TNFRSF1A', 'TNF_TNFRSF1B',
                                                             'TNFSF10_TNFRSF10B', 'TNFSF14_TNFRSF14', 
                                                             'IL34_CSF1R'))),
                 remove.isolate = FALSE, color.heatmap=c('viridis'), n.colors = 3,  line.on = T, grid.on = F,
                 sort.by.source = F, sort.by.target = T, sort.by.source.priority = FALSE) + coord_flip()

ggsave(filename = file.path(path, 'interactions_4.pdf'),
       scale = 0.5, width = 28, height = 16, units='cm')


