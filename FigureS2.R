# ==============================================================================
# Supplementary Figure 2
#   A  Average expression of macrophage marker genes per macrophage state,
#      z-scored per gene, with hierarchical clustering of states and genes
#
# Input from the Figure 2 script:
#   multiome_mac
# ==============================================================================


# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(pheatmap)
  library(RColorBrewer)
})


# ------------------------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------------------------
outdir <- path.expand("~/Desktop/out/figS2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Italic gene labels for pheatmap
italic_labels <- function(x) parse(text = paste0("italic('", x, "')"))

mac_levels <- c("Phago Mac", "Infl Mac", "Endo Mac", "Res Mac")


# ==============================================================================
# Panel A - Heatmap of macrophage markers
# ==============================================================================
# --- Average expression per macrophage state ---
avg <- as.matrix(AverageExpression(multiome_mac, assays = "SCT", group.by = "MacState")$SCT)

# Seurat can alter spaces / underscores in group names, so match them loosely
name_key <- function(x) gsub("[^A-Za-z]", "", x)
idx <- match(name_key(mac_levels), name_key(colnames(avg)))
stopifnot("macrophage states not found in AverageExpression output" = !anyNA(idx))
avg <- avg[, idx]
colnames(avg) <- gsub(" ", "-", mac_levels)  # the figure writes e.g. "Infl-Mac"

# --- Marker genes ---
marker_genes <- c(
  "CD36", "LYVE1", "RNASE1", "NFKBIZ", "CCL8", "F13A1", "CCL2", "CD80",
  "MAF", "CD209", "MRC1", "CXCL10", "CCL4", "CD274", "CXCL9", "IRF1", "NFKB1", "CCL21",
  "ICAM1", "CCL3", "CD44", "CD58", "VCAM1", "IL1B", "CD83", "IRF8", "APOE", "SELENOP",
  "CD163", "CD38", "MMP9", "CXCL12", "PLAU", "C1QA", "MERTK", "IL6"
)
missing_genes <- setdiff(marker_genes, rownames(avg))
if (length(missing_genes) > 0) {
  message("Not in the SCT assay, omitted from the heatmap: ", paste(missing_genes, collapse = ", "))
}
heat_mat <- t(avg[intersect(marker_genes, rownames(avg)), , drop = FALSE])  # states x genes

# --- Plot (z-score per gene across states, then cluster both axes) ---
pheatmap(
  heat_mat,
  scale             = "column",
  color             = brewer.pal(9, "BuPu"),
  clustering_method = "ward.D2",
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = 12,
  fontsize_col      = 12,
  labels_col        = italic_labels(colnames(heat_mat)),
  filename          = file.path(outdir, "figS2A_macrophage_marker_heatmap.pdf"),
  width             = 12,
  height            = 4.6
)
