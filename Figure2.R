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

# Figure 2a - Macrophage diffusion map
# Extract macrophages and pre-process
multiome_mac <- subset(multiome_myeloid, subset=CellType=='Macrophage')
multiome_mac <- subset(multiome_mac, subset=Cluster%in%c('6'), invert=T)

multiome_mac <- CreateSeuratObject(counts = multiome_mac@assays$RNA@counts, meta.data = multiome_mac@meta.data)
multiome_mac <- SCTransform(multiome_mac, vst.flavor = "v2")
multiome_mac <- RunPCA(multiome_mac, npcs = 100)
multiome_mac <- RunHarmony(multiome_mac, group.by.vars='Sex')
multiome_mac <- RunUMAP(multiome_mac, reduction = 'harmony', dims = c(1:10))

pca <- as.data.frame(t(multiome_mac@reductions[["harmony"]]@cell.embeddings))
colnames(pca) <- multiome_mac$Cluster

# Diffuaion map
dm <- DiffusionMap(t(pca[1:20,]), k=30, sigma='local', n_local=10, rotate=F)

# Add to Seurat object
dr <- multiome_mac@reductions[["umap"]]
dr@cell.embeddings[,1] <- as.numeric(dm$DC1)*100
dr@cell.embeddings[,2] <- as.numeric(dm$DC2)*100
dr@cell.embeddings <- cbind(dr@cell.embeddings, as.numeric(dm$DC3)*100)
colnames(dr@cell.embeddings) <- c('DC_1', 'DC_2', 'DC_3')
dr@key <- 'DC_'
multiome_mac@reductions[["DM"]] <- dr

# Cluster
multiome_mac <- FindNeighbors(multiome_mac, reduction = "harmony", dims = 1:20) %>% FindClusters(algorithm=1, resolution = 1.5)

Idents(multiome_mac) <- paste('c', Idents(multiome_mac), sep='')
multiome_mac <- RenameIdents(object = multiome_mac, `c3` = "D")
multiome_mac <- RenameIdents(object = multiome_mac, `c10` = "D")
multiome_mac <- RenameIdents(object = multiome_mac, `c5` = "D")
multiome_mac <- RenameIdents(object = multiome_mac, `c8` = "D")
multiome_mac <- RenameIdents(object = multiome_mac, `c0` = "D")
multiome_mac <- RenameIdents(object = multiome_mac, `c4` = "A")
multiome_mac <- RenameIdents(object = multiome_mac, `c6` = "A")
multiome_mac <- RenameIdents(object = multiome_mac, `c1` = "C")
multiome_mac <- RenameIdents(object = multiome_mac, `c9` = "C")
multiome_mac <- RenameIdents(object = multiome_mac, `c2` = "C")
multiome_mac <- RenameIdents(object = multiome_mac, `c7` = "B")
multiome_mac$Cluster <- Idents(multiome_mac)

# Plot
Idents(multiome_mac) <- factor(Idents(multiome_mac), levels=(c('A', 'B', 'C', 'D')))
DimPlot(multiome_mac, reduction='DM', group.by = 'Cluster', dims = c(1, 2), label = F) + 
  ggtitle('') + NoAxes() + theme_void()

ggsave(filename = file.path(path, 'umap_macrophage.pdf'),
       scale = 0.5, width = 25, height = 20, units='cm')


# Figure 2b - Macrophage proportions
# Format data
multiome_mac_control <- subset(multiome_mac, subset=Condition=='Control')
multiome_mac_uuo <- subset(multiome_mac, subset=Condition=='UUO')

control_df <- data.frame(table(multiome_mac_control$Cluster))
control_df$pct <- control_df$Freq/sum(control_df$Freq)*100
control_df$Var1 <- factor(control_df$Var1, levels = c('A', 'B', 'C', 'D'))
uuo_df <- data.frame(table(multiome_mac_uuo$Cluster))
uuo_df$pct <- uuo_df$Freq/sum(uuo_df$Freq)*100
uuo_df$Var1 <- factor(uuo_df$Var1, levels = c('A', 'B', 'C', 'D'))

# Plot
ggplot(control_df, aes(x="", y=pct, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=c('A'='#00BFC4', 'B'='#F8766D', 'C'='#7CAE00', 'D'='#C77CFF'))

ggsave(filename = file.path(path, 'pie_control.pdf'),
       scale = 0.5, width = 20, height = 15, units='cm')

ggplot(uuo_df, aes(x="", y=pct, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + scale_fill_manual(values=c('A'='#00BFC4', 'B'='#F8766D', 'C'='#7CAE00', 'D'='#C77CFF'))

ggsave(filename = file.path(path, 'pie_uuo.pdf'),
       scale = 0.5, width = 20, height = 15, units='cm')


# Figure 2c - Macrophage pathway scores
go_terms <- readRDS(file.path(path, 'go_terms.rds'))

# Add scores for terms of interest
pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0019221'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(cytokine=pw_genes$SYMBOL)) 
pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0034341'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(IFNG=pw_genes$SYMBOL)) 
multiome_mac@meta.data['response to type II interferon'] <- multiome_mac$IFNG_UCell

pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0034612'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(TNF=pw_genes$SYMBOL)) 
multiome_mac@meta.data['response to tumor necrosis factor'] <- multiome_mac$TNF_UCell

pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0006898'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(endocytosis=pw_genes$SYMBOL)) 
multiome_mac@meta.data['receptor-mediated endocytosis'] <- multiome_mac$endocytosis_UCell

pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0006910'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(phago=pw_genes$SYMBOL)) 
multiome_mac@meta.data['phagocytosis, recognition'] <- multiome_mac$phago_UCell

pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0090130'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(migration=pw_genes$SYMBOL)) 
multiome_mac@meta.data['tissue migration'] <- multiome_mac$migration_UCell

pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0021700'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(development=pw_genes$SYMBOL)) 
multiome_mac@meta.data['developmental maturation'] <- multiome_mac$development_UCell

pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0019882'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(antigen=pw_genes$SYMBOL)) 
multiome_mac@meta.data['antigen processing and presentation'] <- multiome_mac$antigen_UCell

pw_genes <- genes <- bitr(unname(unlist(go_terms['GO:0019915'])), fromType = "ENTREZID", toType = "SYMBOL", OrgDb='org.Hs.eg.db')
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list(lipid=pw_genes$SYMBOL)) 

# Plot
DotPlot_scCustom(seurat_object = multiome_mac, 
                 features=rev(c('phagocytosis, recognition', 'receptor-mediated endocytosis', 'antigen processing and presentation', 'tissue migration',
                                'developmental maturation', 'response to tumor necrosis factor', 'response to type II interferon')), 
                 flip_axes = T, scale.min = 0, scale.max=50) + scale_color_viridis() + NoLegend()

ggsave(filename = file.path(path, 'macrophage_pathways.pdf'),
       scale = 0.5, width = 20, height = 20, units='cm')


# Figure 2d - Macrophage genes dotplot
DotPlot(multiome_mac, 
        features = rev(c('CCL2', 'CCL3', 'CCL4', 'CCL8', 'CD274', 'CXCL9', 'CXCL10', 'IL1B', 'TNFRSF4')), 
        scale.min = 0, 
        cols = c('grey90', 'navy')) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) + NoLegend()

ggsave(filename = file.path(path, 'macrophage_dotplot.pdf'),
       scale = 0.5, width = 11, height = 20, units='cm')


# Figure 2f - Macrophage dotplot of enriched GO terms
counts <- read.csv(file.path(path, "mac_bulkRNA_counts.csv"))
rownames(counts) <- counts$X; counts$X <- NULL
colnames(counts) <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b", "5a", "5b")
sample_sheet <- data.frame(names = c("1a.tab", "1b.tab", "2a.tab", "2b.tab", 
                                            "3a.tab", "3b.tab", "4a.tab", "4b.tab", "5a.tab", "5b.tab"), 
                                  donor = c("M1_1", "M1_1", "M1_2", "M1_2", "M1_3", "M1_3", 
                                            "M1_4", "M1_4", "M1_5", "M1_5"), 
                                  celltype = c("Macrophage", "Macrophage", "Macrophage", "Macrophage", "Macrophage", "Macrophage", 
                                               "Macrophage", "Macrophage", "Macrophage", "Macrophage"), 
                                  condition = c("NIR_RPTEC", "IR_RPTEC", "NIR_RPTEC", "IR_RPTEC", 
                                                "NIR_RPTEC", "IR_RPTEC", "NIR_RPTEC", "IR_RPTEC", "NIR_RPTEC", 
                                                "IR_RPTEC"))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_sheet,
                              design = ~ condition+donor)
dds <- dds[rowSums(counts(dds)) >= 10,]
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

# Process DE results
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC <- as.data.frame(resLFC)
resLFC$gene_id <- sub("\\..*", "", rownames(resLFC))

# Translate gene IDs
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), filters = "ensembl_gene_id", 
                 values = resLFC$gene_id, mart = ensembl) 
resLFC <- merge(resLFC, genemap, by.x="gene_id", by.y="ensembl_gene_id")

# Subset by LFC
resLFC_subset <- resLFC[resLFC$log2FoldChange<(-0.5) & resLFC$padj<0.05,]

# Prepare query genes
geneList = resLFC_subset$log2FoldChange
names(geneList) = as.character(resLFC_subset$entrezgene_id)
geneList = sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(names(geneList))]
geneList <- geneList[!duplicated(names(geneList))]
universe <- as.character(resLFC$entrezgene_id)

ego <- enrichGO(gene          = names(geneList),
                universe      = universe,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.2,
                qvalueCutoff  = 0.2,
                readable      = TRUE)


results <- ego@result
results_subset <- results[results$ID%in%c('GO:0001819', 'GO:0019221', 'GO:0034341', 'GO:0050900', 'GO:0050727', 'GO:0045088', 'GO:0071356', 'GO:0022409'),]
results_subset <- results_subset[order(results_subset$FoldEnrichment, decreasing = T),]
results_subset$GeneRatio_numeric <- sapply(strsplit(results_subset$GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

max_count <- max(results_subset$Count, na.rm = TRUE)

# Plot
ggplot(results_subset, aes(x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_viridis(option = "D", direction = -1, name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count", limits = c(0, max_count)) +
  labs(x = "Gene Ratio", y = "GO Term", title = "GO Term Enrichment") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

ggsave(filename = file.path(path, 'macrophage_go_terms.pdf'),
       scale = 0.5, width = 34, height = 20, units='cm')


# Figure 2f - Macrophage dotplot of enriched GO terms
# Extract count matrix
matrix <- as.data.frame(assay(vsd))
matrix$gene_id <- sub("\\..*", "", rownames(matrix))

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", 
                 values = matrix$gene_id, mart = ensembl) 
matrix <- merge(matrix, genemap, by.x="gene_id", by.y="ensembl_gene_id")
rownames(matrix) <- make.unique(matrix$hgnc_symbol)
matrix$hgnc_symbol <- NULL
matrix$gene_id <- NULL
colnames(matrix) <- sub('\\.', '-', colnames(matrix))

# Prepare plot
genes <- c('MERTK', 'CABLES1', 'CDKN1C', 'PDK4', 'ADAMTS15', 'MS4A6A', 'CEACAM4', 'THBS1', 'THBS1-IT1',
           'CCL1', 'CCL2', 'CCL3', 'CCL3L3', 'CCL4', 'CCL8', 'CD274', 'CXCL9', 'CXCL10', 
           'IL1B', 'IL6', 'IL12B', 'IL17A', 'IL17F', 'TNF', 'TNFSF15', 'TNFRSF4',
           'TGFB1', 'CSF2', 'CLEC5A', 'EDN1', 'HSD11B1', 'SERPINE1', 'SIMALR')

select <- is.element(rownames(matrix), genes)
matrix_subset <- matrix[select,]
matrix_subset <- matrix_subset[genes, c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10), drop = FALSE]

col_data <- data.frame(Condition=c(rep('NIR_RPTEC', 5), rep('IR_RPTEC', 5)))
rownames(col_data) <- c("1a", "2a", "3a", "4a", "5a", "1b", "2b", "3b", "4b", "5b")
cols <- brewer.pal(9, "BuPu")
col_list <- brewer.pal(10, "Paired")
annotation_colors = list(Condition = c('IR_RPTEC'='navy', 'NIR_RPTEC'='grey80'))

# Plot
pheatmap(matrix_subset, cluster_rows=F, cluster_cols=F, show_rownames=T, show_colnames=T, 
         color=cols, 
         clustering_method="ward.D2",
         scale="row", 
         annotation_col=col_data, 
         fontsize_row=12, 
         annotation_colors = annotation_colors,
         labels_row = parse(text = paste0("italic('", rownames(matrix_subset), "')")),
         gaps_col = 5, gaps_row = 9)


# Figure 2g - Projection of bulk signature onto sn data
# Score signature
resLFC_subset <- resLFC[resLFC$log2FoldChange<(-0.5) & resLFC$padj<0.05,]
multiome_mac <- AddModuleScore_UCell(multiome_mac, features = list('upregulated'= resLFC_subset$hgnc_symbol),
                                    chunk.size = 8000, ncores = 10, name='', maxRank=3000)

# Plot
VlnPlot(multiome_mac, features = "upregulated", pt.size = 0) +
  scale_fill_manual(values = c("A" = "#00BFC4", "B" = "#F8766D", "C" = "#7CAE00", "D" = "#C77CFF")) +
  NoLegend()

ggsave(filename = file.path(path, 'signature_violin.pdf'),
       scale = 0.5, width = 13, height = 20, units='cm')



# Figure 2l - Macrophage MERTK expression
Idents(multiome_mac) <- factor(Idents(multiome_mac), levels=rev((c('A', 'B', 'C', 'D'))))
DotPlot(multiome_mac, 
        features = 'MERTK', 
        scale.min = 0, 
        cols = c('grey90', 'navy')) + 
  NoLegend() + 
  RotatedAxis() +
  scale_x_discrete(labels = function(x) parse(text = paste0("italic('", x, "')")))

ggsave(filename = file.path(path, 'dotplot_mertk.pdf'),
       scale = 0.5, width = 8, height = 14, units='cm')

