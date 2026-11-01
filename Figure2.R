# ==============================================================================
# Figure 2
#   A  Macrophage UMAP, coloured by polarisation state
#   B  State proportions in normal vs UUO kidney (pie charts)
#   C  GO-term pathway scores per state
#   D  Marker-gene dot plot per state
#   F  GO enrichment of genes induced by SC CM (bulk RNA-seq of macrophages)
#   G  Heatmap of selected genes (bulk RNA-seq)
#   H  Bulk-derived signature projected onto the single-nucleus macrophages
# ==============================================================================


# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(UCell)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})

set.seed(1)

# ------------------------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------------------------
data_dir             <- #directory with supplementary data
path_myeloid_rds     <- file.path(data_dir, "myeloid_data.rds")
path_go_terms_rds    <- file.path(data_dir, "go_terms.rds")
path_mac_bulk_counts <- file.path(data_dir, "mac_bulkRNA_counts.csv")

outdir <- path.expand("~/Desktop/out/fig2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_cores <- 10  # UCell parallelism

# Plot helper
save_plot <- function(plot, file, width, height) {
  ggsave(file.path(outdir, file), plot = plot,
         width = width, height = height, units = "cm")
}

# Helper: italic gene labels for ggplot axes and pheatmap
italic_labels <- function(x) parse(text = paste0("italic('", x, "')"))


# ------------------------------------------------------------------------------
# Colours
# ------------------------------------------------------------------------------
# Convert an HSV colour (3x1 matrix: h, s, v, each in 0-1) to RGB (0-1)
hsv2rgb <- function(x) {
  h <- x[1, 1]
  s <- x[2, 1]
  v <- x[3, 1]
  C <- s * v
  hdash <- h * 6
  X <- C * (1 - abs(hdash %% 2 - 1))
  if (0 <= hdash & hdash <= 1) RGB1 <- c(C, X, 0)
  if (1 <= hdash & hdash <= 2) RGB1 <- c(X, C, 0)
  if (2 <= hdash & hdash <= 3) RGB1 <- c(0, C, X)
  if (3 <= hdash & hdash <= 4) RGB1 <- c(0, X, C)
  if (4 <= hdash & hdash <= 5) RGB1 <- c(X, 0, C)
  if (5 <= hdash & hdash <= 6) RGB1 <- c(C, 0, X)
  RGB1 + (v - C)
}

# Desaturate a colour: multiply its saturation by p (p = 1 leaves it unchanged,
# smaller p gives a paler colour)
pastellize <- function(x, p) {
  if (is.character(x)) x <- col2rgb(x) / 255
  if (is.numeric(x))   x <- matrix(x, nrow = 3)
  col <- rgb2hsv(x, maxColorValue = 1)
  col[2, 1] <- col[2, 1] * p
  col <- hsv2rgb(col)
  rgb(col[1], col[2], col[3])
}

# Material colour ramps (the same ones used for Figure 1)
oranges <- pal_material("orange", alpha = 1)(10)
reds    <- pal_material("red",    alpha = 1)(10)
indigos <- pal_material("indigo", alpha = 1)(10)
greens  <- pal_material("green",  alpha = 1)(10)

mac_levels <- c("Phago Mac", "Infl Mac", "Endo Mac", "Res Mac")
mac_colors <- c(
  "Phago Mac" = pastellize(oranges[8], 0.5),
  "Infl Mac"  = pastellize(reds[8],    0.8),
  "Endo Mac"  = pastellize(indigos[6], 0.8),
  "Res Mac"   = pastellize(greens[7],  0.8)
)


# ------------------------------------------------------------------------------
# Gene ID mapping (offline, via org.Hs.eg.db)
# ------------------------------------------------------------------------------
annotate_ensembl <- function(ids) {
  ids <- sub("\\..*", "", ids)  # strip Ensembl version suffix
  data.frame(
    gene_id = ids,
    symbol  = unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, keytype = "ENSEMBL",
                                           column = "SYMBOL",   multiVals = "first")),
    entrez  = unname(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, keytype = "ENSEMBL",
                                           column = "ENTREZID", multiVals = "first")),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# Panel A - Macrophage UMAP
# ==============================================================================

# --- Extract macrophages, drop low-quality cluster 6 (from the myeloid object) ---
multiome_myeloid <- readRDS(path_myeloid_rds)
multiome_mac <- subset(multiome_myeloid,
                       subset = CellType == "Macrophage" & !(Cluster %in% "6"))

# --- Rebuild as a fresh object, then normalise, integrate and embed ---
multiome_mac <- CreateSeuratObject(
  counts    = GetAssayData(multiome_mac, assay = "RNA", slot = "counts"),
  meta.data = multiome_mac@meta.data
)
multiome_mac <- SCTransform(multiome_mac, vst.flavor = "v2")
multiome_mac <- RunPCA(multiome_mac, npcs = 100)
multiome_mac <- RunHarmony(multiome_mac, group.by.vars = "Sex")
multiome_mac <- RunUMAP(multiome_mac, reduction = "harmony", dims = 1:30,
                        min.dist = 0.03, n.neighbors = 20)

# --- Cluster ---
multiome_mac <- FindNeighbors(multiome_mac, reduction = "harmony", dims = 1:30)
multiome_mac <- FindClusters(multiome_mac, algorithm = 1, resolution = 1.3)

# --- Cluster markers, used to inform the manual annotation below ---
markers <- SeuratWrappers::RunPrestoAll(multiome_mac, only.pos = TRUE)
write.csv(markers, file.path(outdir, "suppl_macrophage_cluster_markers.csv"), row.names = FALSE)

# --- Manual annotation of clusters into four macrophage states ---
# Specific to this clustering (resolution 1.3, seed 1); re-check against the
# marker table if anything upstream changes.
cluster_to_state <- c(
  "0" = "Phago Mac", "1" = "Endo Mac", "2" = "Res Mac", "3" = "Endo Mac",
  "4" = "Phago Mac", "5" = "Endo Mac", "6" = "Res Mac", "7" = "Infl Mac"
)
stopifnot("unannotated clusters present" =
            all(levels(multiome_mac$seurat_clusters) %in% names(cluster_to_state)))

multiome_mac$MacState <- factor(
  unname(cluster_to_state[as.character(multiome_mac$seurat_clusters)]),
  levels = mac_levels
)
Idents(multiome_mac) <- "MacState"

# --- Plot ---
p_umap <- DimPlot(multiome_mac, group.by = "MacState", cols = mac_colors) +
  ggtitle("") + theme_void()
save_plot(p_umap, "fig2A_umap_macrophages.pdf", width = 9, height = 7)


# ==============================================================================
# Panel B - Macrophage state proportions (Control vs UUO)
# ==============================================================================

state_proportions <- function(obj) {
  df <- as.data.frame(table(obj$MacState))
  df$pct <- df$Freq / sum(df$Freq) * 100
  df
}

plot_pie <- function(df, title) {
  ggplot(df, aes(x = "", y = pct, fill = Var1)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = mac_colors, name = NULL) +
    ggtitle(title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5))
}

p_pie_control <- plot_pie(state_proportions(subset(multiome_mac, subset = Condition == "Control")),
                          "Normal kidney")
p_pie_uuo     <- plot_pie(state_proportions(subset(multiome_mac, subset = Condition == "UUO")),
                          "UUO kidney")

p_pies <- (p_pie_control / p_pie_uuo) + plot_layout(guides = "collect")  # one shared legend
save_plot(p_pies, "fig2B_pies_state_proportions.pdf", width = 11, height = 12)


# ==============================================================================
# Panel C - Macrophage pathway scores (GO terms)
# ==============================================================================

# GO term -> label shown on the plot, in top-to-bottom plotting order
go_sets <- c(
  "GO:0006910" = "phagocytosis, recognition",
  "GO:0006898" = "receptor-mediated endocytosis",
  "GO:0019882" = "antigen processing and presentation",
  "GO:0090130" = "tissue migration",
  "GO:0021700" = "developmental maturation",
  "GO:0034612" = "response to tumor necrosis factor",
  "GO:0034341" = "response to type II interferon"
)

go_terms <- readRDS(path_go_terms_rds)  # list: GO ID -> Entrez IDs

go_symbols <- lapply(names(go_sets), function(go_id) {
  entrez  <- as.character(unname(unlist(go_terms[go_id])))
  symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID",
                                   column = "SYMBOL", multiVals = "first")
  unique(stats::na.omit(unname(symbols)))
})
names(go_symbols) <- sub(":", "_", names(go_sets))  # safe feature names for UCell
stopifnot(all(lengths(go_symbols) > 0))

multiome_mac <- AddModuleScore_UCell(multiome_mac, features = go_symbols,
                                     chunk.size = 8000, ncores = n_cores, name = "")
for (i in seq_along(go_sets)) {  # copy each score to a readable metadata column
  multiome_mac[[go_sets[[i]]]] <- multiome_mac[[names(go_symbols)[i]]][, 1]
}

p_pathways <- DotPlot(multiome_mac, group.by = "MacState", features = rev(unname(go_sets)),
                      scale.min = 0, scale.max = 50) +
  coord_flip() +
  scale_color_viridis_c() +
  guides(size = "none") +
  labs(color = "Score", x = NULL, y = NULL) +
  RotatedAxis()
save_plot(p_pathways, "fig2C_macrophage_pathways.pdf", width = 13, height = 9)


# ==============================================================================
# Panel D - Macrophage marker gene dot plot
# ==============================================================================

p_markers <- DotPlot(
  multiome_mac, group.by = "MacState",
  features = rev(c("CCL2", "CCL3", "CCL4", "CCL8", "CD274", "CXCL9", "CXCL10", "IL1B", "TNFRSF4")),
  scale.min = 0, cols = c("grey90", "navy")
) +
  coord_flip() +
  scale_x_discrete(labels = italic_labels) +
  labs(x = NULL, y = NULL) +
  theme(legend.title = element_blank()) +
  RotatedAxis()
save_plot(p_markers, "fig2D_macrophage_markers.pdf", width = 8, height = 9)


# ==============================================================================
# Panel F - GO enrichment (bulk RNA-seq, macrophages +/- SC CM)
# ==============================================================================

# --- Load counts (first column = gene ID) ---
counts <- read.csv(path_mac_bulk_counts, check.names = FALSE)
rownames(counts) <- counts[[1]]
counts <- counts[, -1, drop = FALSE]
stopifnot("non-numeric columns in count matrix" = all(vapply(counts, is.numeric, logical(1))),
          "expected 10 samples"                 = ncol(counts) == 10)

# Assumes the file's columns are in this order (a = control, b = SC CM, per donor)
sample_ids <- c("1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b", "5a", "5b")
colnames(counts) <- sample_ids

sample_sheet <- data.frame(
  donor     = rep(paste0("M1_", 1:5), each = 2),
  celltype  = "Macrophage",
  condition = rep(c("NIR_RPTEC", "IR_RPTEC"), 5),  # NIR = control, IR = SC CM
  row.names = sample_ids
)

# --- DESeq2 ---
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_sheet,
                              design = ~ condition + donor)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# Condition levels sort alphabetically (IR_RPTEC, NIR_RPTEC), so coefficient 2 is
# NIR vs IR. The sign is flipped below so that positive = higher with SC CM.
stopifnot(resultsNames(dds)[2] == "condition_NIR_RPTEC_vs_IR_RPTEC")
de <- as.data.frame(lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm"))
de <- cbind(annotate_ensembl(rownames(de)), de)
de$log2FoldChange <- -de$log2FoldChange
write.csv(de, file.path(outdir, "suppl_macrophage_bulk_DE.csv"), row.names = FALSE)

# --- Genes induced by SC CM (log2FC > 0.5, padj < 0.05) ---
de_up <- de[!is.na(de$padj) & de$padj < 0.05 & de$log2FoldChange > 0.5, ]
universe_entrez <- unique(stats::na.omit(de$entrez))

ego <- enrichGO(
  gene = unique(stats::na.omit(de_up$entrez)), universe = universe_entrez,
  OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
  pvalueCutoff = 0.2, qvalueCutoff = 0.2, readable = TRUE
)
write.csv(ego@result, file.path(outdir, "suppl_macrophage_go_enrichment.csv"), row.names = FALSE)

# --- GO terms shown in the figure ---
go_plot_ids <- c("GO:0001819", "GO:0019221", "GO:0034341", "GO:0050900",
                 "GO:0050727", "GO:0045088", "GO:0071356", "GO:0022409")
missing_go <- setdiff(go_plot_ids, ego@result$ID)
if (length(missing_go) > 0) stop("GO terms not enriched: ", paste(missing_go, collapse = ", "))

go_plot <- ego@result[ego@result$ID %in% go_plot_ids, ]
go_plot$GeneRatio_numeric <- vapply(strsplit(go_plot$GeneRatio, "/"),
                                    function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))

p_go <- ggplot(go_plot, aes(x = GeneRatio_numeric, y = reorder(Description, GeneRatio_numeric))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_viridis_c(option = "D", direction = -1, name = "Adjusted p-value") +
  scale_size_continuous(name = "Gene Count", limits = c(0, max(go_plot$Count))) +
  labs(x = "Gene Ratio", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y  = element_text(size = 10, color = "black"),
    axis.text.x  = element_text(color = "black"),
    axis.title   = element_text(color = "black"),
    legend.text  = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )
save_plot(p_go, "fig2F_macrophage_go_terms.pdf", width = 17, height = 10)


# ==============================================================================
# Panel G - Heatmap of selected genes (bulk RNA-seq)
# ==============================================================================

# --- VST matrix with gene symbols as row names ---
vsd_mat <- assay(vsd)
vsd_ann <- annotate_ensembl(rownames(vsd_mat))
keep    <- !is.na(vsd_ann$symbol)
vsd_mat <- vsd_mat[keep, ]
rownames(vsd_mat) <- make.unique(vsd_ann$symbol[keep])

# --- Genes of interest: first block = homeostatic / resident, second = inflammatory ---
heatmap_genes <- c(
  "MERTK", "CD163", "CABLES1", "CDKN1C", "PDK4", "ADAMTS15", "MS4A6A", "CEACAM4", "THBS1", "THBS1-IT1",
  "CD80", "VCAM1", "ICAM1", "IRF8", "CD83", "CD58", "CD44", "CCL1", "CCL2", "CCL3", "CCL3L3", "CCL4", "CCL8",
  "CD274", "CXCL9", "CXCL10", "IL1B", "IL6", "IL12B", "IL17A", "IL17F", "TNF", "TNFSF15", "TNFRSF4",
  "TGFB1", "CSF2", "CLEC5A", "EDN1", "HSD11B1", "SERPINE1", "SIMALR"
)
missing_genes <- setdiff(heatmap_genes, rownames(vsd_mat))
if (length(missing_genes) > 0) {
  stop("Genes not found after ID mapping: ", paste(missing_genes, collapse = ", "),
       ". Likely a symbol difference between org.Hs.eg.db and Ensembl (e.g. lncRNAs).")
}

col_order   <- c(paste0(1:5, "a"), paste0(1:5, "b"))  # controls first, then SC CM
heat_mat    <- vsd_mat[heatmap_genes, col_order]

annot_col <- data.frame(Condition = rep(c("Control", "SC CM"), each = 5), row.names = col_order)
annot_colors <- list(Condition = c("SC CM" = "navy", "Control" = "grey80"))

pheatmap(
  heat_mat,
  scale             = "row",
  color             = brewer.pal(9, "BuPu"),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  annotation_col    = annot_col,
  annotation_colors = annot_colors,
  annotation_legend = FALSE,
  fontsize_row      = 12,
  labels_row        = italic_labels(rownames(heat_mat)),
  gaps_col          = 5,
  gaps_row          = 10,
  filename          = file.path(outdir, "fig2G_heatmap_bulk.pdf"),
  width             = 5,
  height            = 10
)


# ==============================================================================
# Panel H - Bulk-derived signature projected onto single-nucleus macrophages
# ==============================================================================

sc_cm_signature <- unique(stats::na.omit(de_up$symbol))
multiome_mac <- AddModuleScore_UCell(
  multiome_mac, features = list(SC_CM_signature = sc_cm_signature),
  chunk.size = 8000, ncores = n_cores, name = "", maxRank = 3000
)

p_vln <- VlnPlot(multiome_mac, features = "SC_CM_signature", group.by = "MacState", pt.size = 0) +
  scale_fill_manual(values = mac_colors) +
  labs(title = NULL, x = NULL, y = "Signature scores") +
  NoLegend()
save_plot(p_vln, "fig2H_signature_violin.pdf", width = 6.5, height = 10)

