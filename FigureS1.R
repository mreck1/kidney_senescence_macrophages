# ==============================================================================
# Supplementary Figure 1 
#   A  Whole-multiome UMAP, coloured by level-1 annotation
#   B  PT marker genes on the PT diffusion map
#   C  PT diffusion map, samples split by condition
#   D  Volcano plot, irradiated vs control RPTECs
#   E  Myeloid lineage UMAP
#   F  Myeloid marker dot plot
#   G  CellChat: senescent PT -> myeloid ligand-receptor interactions
#
# Inputs from the Figure 1 script:
#  multiome_pt
#  dds
# ==============================================================================


# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(harmony)
  library(DESeq2)
  library(EnhancedVolcano)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(ggsci)
})

set.seed(1)


# ------------------------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------------------------
data_dir             <- #directory with supplementary data
path_multiome_rds    <- # path to multiome .rds object
path_myeloid_saved     <- file.path(data_dir, "myeloid_data.rds")            # myeloid object used by Figure 2
path_myeloid_excluded  <- file.path(data_dir, "myeloid_excluded_cells.txt")  # barcodes removed as artefact (see panel E)

fig1_dir     <- path.expand("~/Desktop/out/figS11")

outdir <- path.expand("~/Desktop/out/figS1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_workers <- 10  # CellChat parallelism

require_file <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    stop("File not found: ", path, if (!is.null(hint)) paste0("\n", hint), call. = FALSE)
  }
  invisible(path)
}

# Plot helper
save_plot <- function(plot, file, width, height, dpi = 300) {
  ggsave(file.path(outdir, file), plot = plot,
         width = width, height = height, units = "cm", dpi = dpi)
}


# ------------------------------------------------------------------------------
# Colours
# ------------------------------------------------------------------------------
# Same helpers as in the Figure 2 script (note: they drop any alpha channel).
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

pastellize <- function(x, p) {
  if (is.character(x)) x <- col2rgb(x) / 255
  if (is.numeric(x))   x <- matrix(x, nrow = 3)
  col <- rgb2hsv(x, maxColorValue = 1)
  col[2, 1] <- col[2, 1] * p
  col <- hsv2rgb(col)
  rgb(col[1], col[2], col[3])
}

# Material colour ramps. pastellize() drops alpha, so the alpha = 0.5 ramps the
# original script used for panel A give exactly the same colours as these.
material  <- function(name) pal_material(name, alpha = 1)(10)
purples     <- material("deep-purple")
blues       <- material("blue")
indigos     <- material("indigo")
reds        <- material("red")
greens      <- material("green")
oranges     <- material("orange")
blue_greys  <- material("blue-grey")
light_blues <- material("light-blue")
greys       <- material("grey")
browns      <- material("brown")
pinks       <- material("pink")

colours_multiome_lvl1 <- c(
  "TAL"          = pastellize(indigos[5],     0.8), "T Cell"    = pastellize(browns[7],      0.5),
  "Myeloid Cell" = pastellize(oranges[8],     0.5), "Endothelia" = pastellize(reds[5],       0.8),
  "PT"           = pastellize(purples[5],     0.8), "IC-B"      = pastellize(greens[3],      0.8),
  "IC-A"         = pastellize(greens[7],      0.8), "CNT"       = pastellize(blues[5],       0.8),
  "B Cell"       = pastellize(pinks[5],       0.3), "DCT"       = pastellize(blues[9],       0.8),
  "PC"           = pastellize(blues[2],       0.8), "Interstitium" = pastellize(light_blues[3], 0.8),
  "PEC"          = pastellize(greys[6],       0.3), "ATL"       = pastellize(blue_greys[2],  0.8),
  "Podocyte"     = pastellize(greys[8],       0.3), "DTL"       = pastellize(blue_greys[4],  0.8)
)

myeloid_levels <- c("CD16 Monocyte", "CD14 Monocyte", "Intermediate", "cDC", "Macrophage")
myeloid_colors <- c(
  "Macrophage"    = purples[4],
  "CD14 Monocyte" = oranges[5],
  "CD16 Monocyte" = oranges[8],
  "Intermediate"  = oranges[2],
  "cDC"           = browns[6]
)


# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
multiome    <- readRDS(path_multiome_rds)


# ==============================================================================
# Panel A - Whole multiome dataset UMAP
# ==============================================================================
p_umap <- DimPlot(
  multiome, group.by = "Annotation.Lvl1", reduction = "umap_wnn",
  cols = colours_multiome_lvl1, pt.size = 0.1, order = FALSE, label = FALSE
) + NoLegend() + NoAxes() + ggtitle("")

p_umap <- LabelClusters(p_umap, id = "Annotation.Lvl1", size = 5, fontface = "bold",
                        color = "black", box = FALSE, repel = FALSE)
save_plot(p_umap, "figS1A_multiome_umap.png", width = 18, height = 15)


# ==============================================================================
# Panel B - PT marker genes on the PT diffusion map
# ==============================================================================
stopifnot("DM" %in% Reductions(multiome_pt))

# TRUE = draw expressing cells on top (for the sparse injury markers)
pt_markers <- c(CUBN = FALSE, MME = FALSE, HAVCR1 = TRUE, VCAM1 = TRUE, CCL2 = TRUE, CDKN1A = TRUE)

marker_plots <- lapply(names(pt_markers), function(gene) {
  FeaturePlot(multiome_pt, features = gene, reduction = "DM",
              order = pt_markers[[gene]], cols = c("grey90", "navy")) +
    NoAxes() +
    theme(plot.title = element_text(face = "italic", hjust = 0.5))
})
p_pt_markers <- wrap_plots(marker_plots, ncol = 3)
save_plot(p_pt_markers, "figS1B_pt_markers.png", width = 37.5, height = 20)


# ==============================================================================
# Panel C - PT diffusion map, samples split by condition
# ==============================================================================
p_dm_samples <- DimPlot(multiome_pt, reduction = "DM", group.by = "Sample",
                        split.by = "Condition", shuffle = TRUE, seed = 1) +
  ggtitle("") + NoAxes()
save_plot(p_dm_samples, "figS1C_dm_samples_by_condition.png", width = 15, height = 10)


# ==============================================================================
# Panel D - RPTEC irradiation volcano plot
# ==============================================================================
# Condition levels sort alphabetically (IR, NR), so coefficient 2 is NR vs IR;
# the sign is flipped so that irradiated > control is positive.
stopifnot(resultsNames(dds)[2] == "Condition_NR_vs_IR")
deg <- as.data.frame(lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm"))
deg$log2FoldChange <- -deg$log2FoldChange
deg <- deg[!is.na(deg$pvalue), ]

# Gene counts shown on the plot (same thresholds as the colouring: raw p < 0.05, |log2FC| > 1)
n_control    <- sum(deg$pvalue < 0.05 & deg$log2FoldChange < -1)
n_irradiated <- sum(deg$pvalue < 0.05 & deg$log2FoldChange >  1)
message("Volcano: ", n_control, " genes higher in control, ",
        n_irradiated, " genes higher in irradiated.")

p_volcano <- EnhancedVolcano(
  deg, lab = rownames(deg), selectLab = "",  # no gene labels
  x = "log2FoldChange", y = "pvalue",
  pCutoff = 0.05, FCcutoff = 1,
  title = "", subtitle = "", caption = "", legendPosition = "none",
  col = c("grey60", "grey60", "grey60", "red4"), colAlpha = 1
) +
  annotate("text", x = -Inf, y = Inf, label = paste(n_control, "genes"),
           hjust = -0.2, vjust = 1.5, size = 4) +
  annotate("text", x = Inf, y = Inf, label = paste(n_irradiated, "genes"),
           hjust = 1.2, vjust = 1.5, size = 4) +
  theme(axis.text = element_text(size = 10, color = "black"))
save_plot(p_volcano, "figS1D_rptec_volcano.png", width = 10, height = 10)


# ==============================================================================
# Panels E and F - Myeloid lineage: UMAP and marker dot plot
# ==============================================================================

# --- Subset myeloid cells (drop pDC and mast cells) and re-process ---
multiome_myeloid <- subset(
  multiome,
  subset = Annotation.Lvl1 %in% "Myeloid Cell" & !(Annotation.Lvl2 %in% c("pDC", "Mast Cell"))
)
multiome_myeloid <- SCTransform(multiome_myeloid, vst.flavor = "v2")
multiome_myeloid <- RunPCA(multiome_myeloid, npcs = 100)
multiome_myeloid <- RunHarmony(multiome_myeloid, assay.use = "SCT", group.by.vars = "Library")
multiome_myeloid <- RunUMAP(multiome_myeloid, reduction = "harmony", dims = 1:50,
                            min.dist = 0.2, n.neighbors = 50)

# --- Cluster ---
multiome_myeloid <- FindNeighbors(multiome_myeloid, reduction = "harmony", dims = 1:50)
multiome_myeloid <- FindClusters(multiome_myeloid, algorithm = 3, resolution = 2)

# --- Annotate ---
# Both maps are keyed on the ORIGINAL FindClusters IDs (c0-c17). The cell-type map
# was previously applied to the already-collapsed cluster labels, whose IDs mean
# something different, so it is applied to the original IDs here.
original_cluster <- paste0("c", multiome_myeloid$seurat_clusters)

# Original cluster -> numbered cluster (18 -> 12), used for macrophage sub-clusters
cluster_map <- c(
  c8 = "7", c13 = "7", c6 = "8", c5 = "8", c10 = "9", c15 = "9",
  c0 = "1", c1 = "2", c3 = "3", c9 = "3", c12 = "4",
  c2 = "5", c17 = "5", c11 = "6", c7 = "10", c4 = "10", c14 = "11", c16 = "12"
)
# Original cluster -> cell type
celltype_map <- c(
  c8 = "CD16 Monocyte", c13 = "CD16 Monocyte", c6 = "CD14 Monocyte",
  c5 = "Intermediate", c10 = "Intermediate", c15 = "Intermediate", c11 = "Intermediate",
  c2 = "Macrophage", c17 = "Macrophage", c0 = "Macrophage", c1 = "Macrophage",
  c3 = "Macrophage", c9 = "Macrophage", c12 = "Macrophage",
  c7 = "cDC", c4 = "cDC", c14 = "cDC", c16 = "cDC"
)
stopifnot("clusters without annotation" =
            all(original_cluster %in% names(cluster_map)) &&
            all(original_cluster %in% names(celltype_map)))

multiome_myeloid$Cluster  <- unname(cluster_map[original_cluster])
multiome_myeloid$CellType <- factor(unname(celltype_map[original_cluster]), levels = myeloid_levels)


# --- Panel E: UMAP ---
p_myeloid <- DimPlot(
  multiome_myeloid, group.by = "CellType", reduction = "umap",
  cols = myeloid_colors, pt.size = 0.8, order = FALSE
) +
  ggtitle("") + NoAxes() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = NA))
save_plot(p_myeloid, "figS1E_myeloid_umap.pdf", width = 15, height = 10)

# --- Panel F: marker dot plot (first level plots at the bottom, so reverse) ---
multiome_myeloid$CellType_dot <- factor(multiome_myeloid$CellType, levels = rev(myeloid_levels))

p_myeloid_dot <- DotPlot(
  multiome_myeloid, group.by = "CellType_dot",
  features = c("FCGR3A", "VCAN", "FCN1", "SORL1", "ITPR1", "FLT3", "IDO1", "MRC1", "CD163"),
  scale.min = 0, cols = c("grey90", "navy")
) +
  RotatedAxis() +
  labs(x = NULL, y = NULL, size = "% Cells") +
  guides(size = guide_legend(override.aes = list(shape = 21, fill = "white"))) +
  theme(axis.text.x = element_text(face = "italic"))
save_plot(p_myeloid_dot, "figS1F_myeloid_dotplot.pdf", width = 15, height = 10)


# ==============================================================================
# Panel G - CellChat: senescent PT -> myeloid signalling
# ==============================================================================

# --- Cell labels: myeloid cell types and PT classes ---
# `class` is a metadata column carried over from the multiome object; it must
# contain a "PT senescent" label.
pt_label_col <- "class"
pt_source    <- "PT senescent"
stopifnot(pt_label_col %in% colnames(multiome_pt@meta.data))

cell_labels <- c(
  setNames(as.character(multiome_myeloid$CellType),         colnames(multiome_myeloid)),
  setNames(as.character(multiome_pt@meta.data[[pt_label_col]]), colnames(multiome_pt))
)
stopifnot("source label not found" = pt_source %in% cell_labels,
          "cells missing from multiome" = all(names(cell_labels) %in% colnames(multiome)))

data <- subset(multiome, cells = names(cell_labels))
data$CellType <- factor(cell_labels[colnames(data)])
stopifnot(!anyNA(data$CellType))

# --- Run CellChat ---
cellchat <- createCellChat(object = GetAssayData(data, assay = "SCT", slot = "data"),
                           meta = data@meta.data, group.by = "CellType")
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = n_workers)
options(future.globals.maxSize = 6000 * 1024^2)

cellchat <- identifyOverExpressedGenes(cellchat, min.cells = 3)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = FALSE,
                              type = "truncatedMean", trim = 0.005)
cellchat <- filterCommunication(cellchat, min.cells = 2)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
future::plan("sequential")

# --- Selected ligand-receptor pairs ---
lr_pairs <- c(
  "CCL2_CCR2", "CCL20_CCR6", "CXCL12_CXCR4", "ICAM1_ITGAX_ITGB2", "VCAM1_INTEGRIN_ADB2",
  "TNF_TNFRSF1A", "TNF_TNFRSF1B", "TNFSF10_TNFRSF10B", "TNFSF14_TNFRSF14", "IL34_CSF1R"
)
missing_lr <- setdiff(lr_pairs, CellChatDB.human$interaction$interaction_name)
if (length(missing_lr) > 0) {
  stop("Not in CellChatDB.human: ", paste(missing_lr, collapse = ", "))
}

# Targets are given by name. The first target plots at the bottom after coord_flip.
p_cellchat <- netVisual_bubble(
  cellchat,
  sources.use = pt_source,
  targets.use = c("CD16 Monocyte", "CD14 Monocyte", "Intermediate", "cDC", "Macrophage"),
  thresh = 0.01,
  pairLR.use = data.frame(interaction_name = lr_pairs),
  remove.isolate = FALSE, color.heatmap = "viridis", n.colors = 3,
  line.on = TRUE, grid.on = FALSE,
  sort.by.source = FALSE, sort.by.target = TRUE, sort.by.source.priority = FALSE
) +
  coord_flip() +
  theme(axis.text.x = element_text(face = "italic"))
save_plot(p_cellchat, "figS1G_cellchat_pt_to_myeloid.pdf", width = 14, height = 8)
