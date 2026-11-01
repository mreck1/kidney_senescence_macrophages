# ==============================================================================
# Figure 1
#   B  Proximal tubule (PT) diffusion map, coloured by PT state
#   C  Diffusion map, coloured by pseudotime
#   F  Heatmap of SASP / cell-cycle genes in irradiated vs control RPTECs
#   G  Senescence-associated signature scores along pseudotime
#   H  Signature-score-based clusters along pseudotime
#   I  Marker dot plots for the Senescent cluster (SASP, signalling, ECM, other)
#   J  CDKN1A chromatin accessibility (ATAC) and expression
# ==============================================================================


# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(harmony)
  library(UCell)
  library(destiny)
  library(DESeq2)
  library(SuperCell)
  library(pheatmap)
  library(ggplot2)
  library(patchwork)
  library(ggsci)
  library(scales)
  library(nichenetr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(1)


# ------------------------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------------------------
data_dir             <- #directory with supplementary data
path_multiome_rds    <- # path to multiome .rds object
path_rptec_counts    <- file.path(data_dir, "rptec_irr_counts.csv")
path_universal_aging    <- file.path(data_dir, "kidney_rodent.csv")
path_tage_results    <- file.path(data_dir, "annotated_tAge_taco_metadata_2026-07-08.csv")

outdir <- path.expand("~/Desktop/out/fig1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

path_taco_expr <- file.path(outdir, "taco_expression_matrix.csv")
path_taco_meta <- file.path(outdir, "taco_metadata.csv")

n_cores       <- 10     # UCell parallelism
dpt_root_cell <- 2746   # index of the DPT root cell within the PT object
mc_gamma      <- 30     # SuperCell graining level (tune to cell counts)
k_clusters    <- 8      # k-means clusters on the combined signature scores

# Plot helper
save_plot <- function(plot, file, width, height) {
  ggsave(file.path(outdir, file), plot = plot,
         width = width, height = height, units = "cm")
}


# ------------------------------------------------------------------------------
# Colour palettes
# ------------------------------------------------------------------------------
indigos <- pal_material("indigo", alpha = 1)(10)
oranges <- pal_material("orange", alpha = 1)(10)

pt_colors <- c(
  "PT S1"             = indigos[5],
  "PT S2"             = indigos[3],
  "PT S3"             = indigos[8],
  "Injury-associated" = oranges[4]
)

cluster_colors <- c(
  "Normal"       = "#A8D9B8",
  "Acute injury" = "#9491C4",
  "Senescent"    = "#D67A7A"
)

score_colors <- c(
  "PTEC Irradiation"   = "#C9852E",
  "tAge"               = "#3D5A80",
  "SenMayo"            = "#7B5EA7",
  "Tabula muris senis" = "#5DA88A"
)

map_ensembl_to_symbol <- function(ids) {
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = ids, keytype = "ENSEMBL",
    column = "SYMBOL", multiVals = "first"
  )
  out <- data.frame(ensembl_gene_id = names(symbols),
                    hgnc_symbol     = unname(symbols))
  out[!is.na(out$hgnc_symbol), ]
}


# ------------------------------------------------------------------------------
# Gene signatures
# ------------------------------------------------------------------------------
sen_mayo <- c(
  "ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6",
  "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3",
  "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2",
  "CSF2RB", "CST4", "CTNNB1", "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16",
  "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG",
  "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF",
  "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4",
  "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18", "IL1A", "IL1B",
  "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA",
  "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14",
  "MMP2", "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF",
  "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4",
  "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2",
  "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA",
  "VEGFC", "VGF", "WNT16", "WNT2"
)

# Tabula Muris Senis signature (mouse symbols; converted to human below)
tabula_muris_mouse <- c(
  "0910001L09Rik", "2010107E04Rik", "2310044H10Rik", "2900097C17Rik",
  "9530068E07Rik", "Aco2", "Actb", "Actg1", "Actr3", "Adam10", "Add1", "Adrm1",
  "Aes", "Aff4", "Amfr", "Anapc11", "Anxa2", "Ap2s1", "Aplp2", "Apoe", "App",
  "Arcn1", "Arhgdia", "Arpc1b", "Arpc3", "Asah1", "Atp1a1", "Atp1b3", "Atp5b",
  "Atp5e", "B2m", "Bri3", "Bsg", "Btg1", "Btg2", "Bzw1", "Cab39", "Calm2",
  "Calr", "Calu", "Canx", "Capn2", "Caprin1", "Capza1", "Capza2", "Cat", "Ccni",
  "Ccpg1", "Cd164", "Cd63", "Cd74", "Cd9", "Cdc42", "Cfl1", "Chmp2a", "Cirbp",
  "Clic1", "Clic4", "Clk1", "Clptm1l", "Cltc", "Cmtm6", "Copa", "Copb1", "Cotl1",
  "Cox6a1", "Cox6b1", "Cox6c", "Cox8a", "Csde1", "Cst3", "Ctnnb1", "Ctsd",
  "Cyba", "Cyfip1", "Dcn", "Ddx3x", "Ddx5", "Ddx6", "Dnaja2", "Dnm2", "Dpysl2",
  "Drap1", "Dusp1", "Eef1a1", "Eef2", "Egr1", "Eif3a", "Eif3b", "Eif3f",
  "Eif4a2", "Eif4g2", "Eif5a", "Erbb2ip", "Errfi1", "Fau", "Fos", "Fth1", "Ftl1",
  "Gadd45gip1", "Galnt1", "Gapdh", "Gas6", "Gdi2", "Ghitm", "Glg1", "Glud1",
  "Glul", "Gm1821", "Gnai3", "Gnb1", "Gnb2l1", "Golph3", "Gpbp1", "Gpx3", "Gpx4",
  "Grn", "Gsk3b", "Gsn", "Gstm1", "H2-Aa", "H2-Ab1", "H2-D1", "H2-K1", "Hbp1",
  "Hif1a", "Hist1h4d", "Hnrnpa2b1", "Hnrnpf", "Hnrnph1", "Hnrnpk", "Hnrnpu",
  "Hsd17b4", "Hspa5", "Hspa8", "Hspa9", "Ier2", "Ifitm2", "Iqgap1", "Itfg1",
  "Itgb1", "Ivns1abp", "Jmjd1c", "Jun", "Junb", "Jund", "Kdm6b", "Klf13",
  "Lamp1", "Lamp2", "Laptm4a", "Lars2", "Lrrc58", "Lrrc8a", "Lyz2", "M6pr",
  "Macf1", "Malat1", "Map1lc3a", "Map3k1", "Mapk1", "Mat2a", "Matr3", "Mbnl1",
  "Mcl1", "Msn", "Mtpn", "Myadm", "Ndufa1", "Nedd4", "Nfic", "Nono", "Npm1",
  "Nptn", "Nucb1", "Nufip2", "Oaz1", "Ogdh", "Osbpl9", "P4hb", "Pan3", "Pbrm1",
  "Pdcd4", "Pdcd6ip", "Pdia3", "Pdia4", "Pebp1", "Pfn1", "Ppia", "Ppp2ca",
  "Ppp2r1a", "Prkar1a", "Pros1", "Prr13", "Prrc2c", "Psap", "Psmd2", "Psmd4",
  "Ptma", "Ptms", "Ptpra", "Pttg1ip", "Rab18", "Rab6a", "Rabac1", "Rac1",
  "Rap1b", "Reep5", "Rhoa", "Rhoc", "Rlim", "Rnaset2b", "Rod1", "Rpl13",
  "Rpl13a", "Rpl14", "Rpl17", "Rpl23", "Rpl23a", "Rpl3", "Rpl31", "Rpl32",
  "Rpl36", "Rpl37", "Rpl37a", "Rpl39", "Rpl4", "Rpl41", "Rpl6", "Rpl9", "Rplp0",
  "Rplp1", "Rpn1", "Rpn2", "Rps10", "Rps12", "Rps14", "Rps16", "Rps17", "Rps18",
  "Rps19", "Rps20", "Rps21", "Rps23", "Rps25", "Rps26", "Rps27", "Rps28",
  "Rps29", "Rps3", "Rps3a", "Rps5", "Rps7", "Rpsa", "Rrp1", "Sdc4", "Sdcbp",
  "Sdha", "Sdhc", "Sec61a1", "Sepp1", "Sept2", "Sepw1", "Serbp1", "Serinc1",
  "Serinc3", "Serp1", "Sf3b2", "Sf3b4", "Sh3bgrl", "Shfm1", "Shisa5", "Skp1a",
  "Slc39a1", "Slc44a2", "Snrpc", "Sod1", "Sparc", "Spnb2", "Spr", "Srsf5",
  "Ssbp4", "Ssr1", "Ssr3", "Stt3a", "Sun2", "Surf4", "Swi5", "Sypl", "Tagln2",
  "Taok1", "Tax1bp1", "Tex261", "Tm9sf2", "Tm9sf3", "Tmbim6", "Tmed9",
  "Tmem123", "Tmem160", "Tmem30a", "Tmsb10", "Tmsb4x", "Tnks2", "Tomm6",
  "Tomm7", "Top2b", "Tpm4", "Tpp1", "Tpt1", "Tram1", "Trf", "Tsc22d3", "Txnip",
  "Uba1", "Uba2", "Uba52", "Ubb", "Ubc", "Ube2m", "Ubxn1", "Usmg5", "Vcp",
  "Vps35", "Wbp2", "Wdr1", "Ybx1", "Ywhaz", "Zfp106", "Zfp207", "Zfp36"
)
tabula_muris <- nichenetr::convert_mouse_to_human_symbols(tabula_muris_mouse)
tabula_muris <- tabula_muris[!is.na(tabula_muris)]


# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
multiome <- readRDS(path_multiome_rds)

# ==============================================================================
# Panels B and C - Proximal tubule diffusion map
# ==============================================================================

# --- Extract PT cells and rebuild as a fresh Seurat object ---
multiome_pt <- subset(multiome, subset = Annotation.Lvl1 %in% "PT")
multiome_pt <- CreateSeuratObject(
  counts    = GetAssayData(multiome_pt, assay = "RNA", slot = "counts"),
  meta.data = multiome_pt@meta.data
)

multiome_pt <- SCTransform(multiome_pt, vst.flavor = "v2")
multiome_pt <- RunPCA(multiome_pt, npcs = 100)
multiome_pt <- RunHarmony(multiome_pt, group.by.vars = "Sample")
DefaultAssay(multiome_pt) <- "SCT"  # explicit: all scoring / DE below uses the SCT assay

# --- Collapse injury / inflammatory PT states into one label for plotting ---
multiome_pt$celltype <- multiome_pt$Annotation.Lvl2
multiome_pt$celltype[multiome_pt$celltype %in% c("PT Injured", "PT Inflammatory")] <- "Injury-associated"

# --- Diffusion map on the first 30 Harmony-corrected PCs ---
harmony_pcs <- Embeddings(multiome_pt, "harmony")[, 1:30]
dm  <- DiffusionMap(harmony_pcs)
dpt <- DPT(dm, tips = dpt_root_cell)
message("DPT root cell: ", colnames(multiome_pt)[dpt_root_cell])

# Store DC1/DC2 as a Seurat reduction (x100 is cosmetic, for plot scaling)
dc_embedding <- cbind(DC_1 = as.numeric(dm$DC1), DC_2 = as.numeric(dm$DC2)) * 100
rownames(dc_embedding) <- colnames(multiome_pt)
multiome_pt[["DM"]] <- CreateDimReducObject(
  embeddings = dc_embedding, key = "DC_", assay = DefaultAssay(multiome_pt)
)

# --- Pseudotime: rank-normalised and flipped so the trajectory runs 0 -> 1 ---
multiome_pt$pseudotime <- 1 - rank(dpt$dpt) / length(dpt$dpt)

# --- Panel B: PT state ---
p_pt <- DimPlot(multiome_pt, reduction = "DM", group.by = "celltype", cols = pt_colors) +
  ggtitle("") + NoAxes()
save_plot(p_pt, "fig1B_dm_pt_states.pdf", width = 12.5, height = 10)

# --- Panel C: pseudotime ---
p_pseudo <- FeaturePlot(multiome_pt, features = "pseudotime", reduction = "DM", order = FALSE) +
  scale_color_viridis_c() + ggtitle("") + NoAxes()
save_plot(p_pseudo, "fig1C_dm_pseudotime.pdf", width = 12.5, height = 10)

# --- Supplementary: sample identity split by condition ---
p_cond <- DimPlot(multiome_pt, reduction = "DM", group.by = "Sample",
                  split.by = "Condition", shuffle = TRUE) +
  ggtitle("") + NoAxes()
save_plot(p_cond, "suppl_dm_sample_by_condition.pdf", width = 15, height = 10)


# ==============================================================================
# Panel F - Heatmap of SASP genes in irradiated RPTECs
# ==============================================================================

# --- Load counts, run DESeq2 ---
counts <- read.csv(path_rptec_counts, row.names = 1)
counts <- counts[rowSums(counts) > 0, ]
rownames(counts) <- make.unique(gsub("\\..*", "", rownames(counts)))  # strip Ensembl version

meta <- data.frame(
  names     = c("B1-NR", "B1-IR", "B3-NR", "B3-IR", "B4-NR", "B4-IR", "B5-NR", "B5-IR", "B6-NR", "B6-IR"),
  Replicate = c("B1", "B1", "B3", "B3", "B4", "B4", "B5", "B5", "B6", "B6"),
  Condition = c("NR", "IR", "NR", "IR", "NR", "IR", "NR", "IR", "NR", "IR")
)
rownames(meta) <- meta$names

# read.csv turns "-" into "." in column names; confirm columns line up with the design table
stopifnot(identical(sub("\\.", "-", colnames(counts)), meta$names))
colnames(counts) <- meta$names

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Condition + Replicate)
dds <- dds[rowSums(counts(dds)) >= 100, ]
dds <- DESeq(dds)
vst_counts <- DESeq2::vst(dds, blind = FALSE)

# --- Map Ensembl IDs to gene symbols ---
rptec_matrix <- as.data.frame(assay(vst_counts))
rptec_matrix$gene_id <- rownames(rptec_matrix)
rptec_matrix <- merge(rptec_matrix, map_ensembl_to_symbol(rownames(rptec_matrix)),
                      by.x = "gene_id", by.y = "ensembl_gene_id")
rownames(rptec_matrix) <- make.unique(rptec_matrix$hgnc_symbol)
rptec_matrix$hgnc_symbol <- NULL
rptec_matrix$gene_id     <- NULL

# --- Subset to genes of interest; order rows and columns (controls, then irradiated) ---
heatmap_genes <- c("CDKN1A", "MKI67", "CCL2", "CCL28", "CXCL1", "CXCL2",
                   "CXCL3", "CXCL6", "CXCL8", "CXCL16")
col_order <- c("B1-NR", "B3-NR", "B4-NR", "B5-NR", "B6-NR",
               "B1-IR", "B3-IR", "B4-IR", "B5-IR", "B6-IR")

stopifnot(all(heatmap_genes %in% rownames(rptec_matrix)))
heat_mat <- rptec_matrix[heatmap_genes, col_order]

annot_col <- data.frame(Treatment = rep(c("Control", "Irradiated"), each = 5),
                        row.names = col_order)
annot_colors <- list(Treatment = c(Irradiated = "purple4", Control = "grey70"))

pheatmap(
  heat_mat,
  scale             = "row",
  color             = colorRampPalette(c(muted("navy", l = 30, c = 70), "white",
                                         muted("red", l = 40, c = 90)))(500),
  annotation_col    = annot_col,
  annotation_colors = annot_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_row          = 2,
  gaps_col          = 5,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize          = 12,
  labels_row        = parse(text = paste0("italic('", rownames(heat_mat), "')")),
  border_color      = "grey50",
  filename          = file.path(outdir, "heatmap_rptec_sasp.pdf"),
  width             = 6,
  height            = 4
)


# ==============================================================================
# Panels G and H - Senescence signatures along pseudotime
# ==============================================================================

# --- 1. Irradiation-response signature from the RPTEC DE results ---
# Condition levels sort alphabetically (IR, NR), so coefficient 2 is NR vs IR;
# the sign is flipped below so that irradiated > control is positive.
stopifnot(resultsNames(dds)[2] == "Condition_NR_vs_IR")
deg <- as.data.frame(lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm"))
deg$gene_id <- rownames(deg)
deg <- merge(deg, map_ensembl_to_symbol(rownames(deg)),
             by.x = "gene_id", by.y = "ensembl_gene_id")
deg$hgnc_symbol     <- make.unique(deg$hgnc_symbol)
deg$log2FoldChange  <- -deg$log2FoldChange
deg <- deg[deg$hgnc_symbol %in% rownames(multiome_pt), ]

genes_up <- deg$hgnc_symbol[which(deg$log2FoldChange > 1 & deg$padj < 0.05)]

# --- 2. Score all three signatures (UCell) ---
score_signature <- function(obj, genes, name) {
  AddModuleScore_UCell(obj, features = setNames(list(genes), name),
                       chunk.size = 8000, ncores = n_cores, name = "")
}
multiome_pt <- score_signature(multiome_pt, genes_up,     "Irradiation")
multiome_pt <- score_signature(multiome_pt, sen_mayo,     "SenMayo")
multiome_pt <- score_signature(multiome_pt, tabula_muris, "TabulaMurisSenis")

# --- 3. Build metacells (SuperCell) ---
multiome_pt <- NormalizeData(multiome_pt)
multiome_pt <- FindVariableFeatures(multiome_pt, nfeatures = 2000)
sc_meta <- multiome_pt@meta.data

SC <- SCimplify(X = GetAssayData(multiome_pt, slot = "data"),
                genes.use = VariableFeatures(multiome_pt),
                gamma = mc_gamma, n.pc = 10)

MC.ge <- supercell_GE(ge = GetAssayData(multiome_pt, slot = "counts"), groups = SC$membership)
MC_celltype <- supercell_assign(cluster = sc_meta$Annotation.Lvl2,
                                supercell_membership = SC$membership, method = "absolute")
MC_sample   <- supercell_assign(cluster = sc_meta$Sample,
                                supercell_membership = SC$membership, method = "absolute")

# --- 4. Average signature scores and pseudotime per metacell ---
score_mat <- t(as.matrix(sc_meta[, c("SenMayo", "TabulaMurisSenis", "Irradiation", "pseudotime")]))
score_mat_mc <- as.data.frame(t(supercell_GE(ge = score_mat, groups = SC$membership, mode = "average")))

# --- 5. Metacell Seurat object with signatures attached ---
MC.seurat <- CreateSeuratObject(counts = MC.ge)
MC.seurat$Annotation.Lvl2  <- as.character(MC_celltype)
MC.seurat$Sample           <- as.character(MC_sample)
MC.seurat$pseudotime       <- score_mat_mc$pseudotime
MC.seurat$Irradiation      <- score_mat_mc$Irradiation
MC.seurat$SenMayo          <- score_mat_mc$SenMayo
MC.seurat$TabulaMurisSenis <- score_mat_mc$TabulaMurisSenis

MC.seurat <- NormalizeData(MC.seurat)
MC.seurat <- FindVariableFeatures(MC.seurat)
MC.seurat <- ScaleData(MC.seurat)
MC.seurat <- RunPCA(MC.seurat, npcs = 80)
MC.seurat <- RunUMAP(MC.seurat, dims = 1:80)
Idents(MC.seurat) <- "Annotation.Lvl2"

# --- 6. Export metacell counts + metadata for TACO ---
expr_mat <- GetAssayData(MC.seurat, slot = "counts")
meta_mc  <- MC.seurat@meta.data
stopifnot(all(colnames(expr_mat) == rownames(meta_mc)))

expr_mat <- expr_mat[rowSums(expr_mat) > 10, ]  # drop genes with negligible counts

write.csv(data.frame(Gene = rownames(expr_mat), as.matrix(expr_mat), check.names = FALSE),
          path_taco_expr, row.names = FALSE)
write.csv(data.frame(SampleID = rownames(meta_mc), meta_mc, check.names = FALSE),
          path_taco_meta, row.names = FALSE)

# --- 7. Import predicted tAge from TACO ---
if (!file.exists(path_tage_results)) {
  stop("TACO results not found at ", path_tage_results, ".\n",
       "Run TACO on:\n  ", path_taco_expr, "\n  ", path_taco_meta,
       "\nthen set `path_tage_results` and re-run.")
}
tage_col <- "tAge.chronological.rodents_kidney.yugene.EN_normalized"
tage <- read.csv(path_tage_results)
if ("SampleID" %in% names(tage)) {
  MC.seurat$tAge <- tage[[tage_col]][match(colnames(MC.seurat), tage$SampleID)]
} else {
  stopifnot(nrow(tage) == ncol(MC.seurat))  # fall back to row order
  MC.seurat$tAge <- tage[[tage_col]]
}

# --- 8. Cluster metacells on the combined ageing / damage signature ---
vars <- c("tAge", "Irradiation", "SenMayo", "TabulaMurisSenis")

set.seed(42)
km <- kmeans(scale(MC.seurat@meta.data[, vars]), centers = k_clusters, nstart = 25)
MC.seurat$score_cluster <- as.character(km$cluster)

# --- 9. Project metacell tAge and cluster back onto single cells ---
# Irradiation / SenMayo / TabulaMurisSenis stay as single-cell UCell scores.
membership <- paste0("Cell_", as.character(SC$membership))
names(membership) <- names(SC$membership)
cell_to_mc <- membership[colnames(multiome_pt)]

mc_idx <- match(cell_to_mc, colnames(MC.seurat))
stopifnot(!anyNA(mc_idx))   # every cell must map to a metacell

multiome_pt$score_cluster <- unname(MC.seurat$score_cluster[mc_idx])
multiome_pt$tAge          <- unname(MC.seurat$tAge[mc_idx])

# --- Panel G: signature scores vs pseudotime ---
plot_df <- multiome_pt@meta.data[, c("pseudotime", vars)]
plot_df[, vars] <- scale(plot_df[, vars])

plot_long <- tidyr::pivot_longer(plot_df, cols = tidyr::all_of(vars),
                                 names_to = "Variable", values_to = "ScaledScore")
plot_long$Variable <- factor(
  plot_long$Variable,
  levels = c("Irradiation", "tAge", "SenMayo", "TabulaMurisSenis"),
  labels = c("PTEC Irradiation", "tAge", "SenMayo", "Tabula muris senis")
)

p_scores <- ggplot(plot_long, aes(x = pseudotime, y = ScaledScore, color = Variable)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.4) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.3, span = 0.3) +
  scale_color_manual(values = score_colors) +
  labs(x = "Pseudotime", y = "Signature score", color = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.title       = element_text(size = 10),
    axis.text        = element_text(color = "grey30"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    legend.position  = "top",
    legend.text      = element_text(size = 11),
    legend.key.width = unit(1.2, "cm"),
    plot.margin      = margin(12, 16, 12, 12)
  )
save_plot(p_scores, "fig1G_scores_vs_pseudotime.pdf", width = 10, height = 5.5)


# --- Panel H: collapse the k-means clusters into three states ---
# Cluster-to-state assignment is specific to set.seed(42) and k = 8.
state_map <- c("4" = "Senescent",    "8" = "Senescent",
               "2" = "Acute injury", "3" = "Acute injury",
               "1" = "Normal", "5" = "Normal", "6" = "Normal", "7" = "Normal")
multiome_pt$senescence_cluster <- unname(state_map[multiome_pt$score_cluster])

plot_df <- data.frame(
  pseudotime         = multiome_pt$pseudotime,
  senescence_cluster = factor(multiome_pt$senescence_cluster,
                              levels = c("Senescent", "Acute injury", "Normal"))  # last level plots on top
)

p_clusters <- ggplot(plot_df, aes(x = pseudotime, y = senescence_cluster, fill = senescence_cluster)) +
  geom_violin(orientation = "y", scale = "width", trim = TRUE,
              color = "black", linewidth = 0.3, alpha = 0.85) +
  geom_jitter(height = 0.15, width = 0, size = 0.3, alpha = 0.15, color = "black") +
  scale_fill_manual(values = cluster_colors) +
  labs(x = "Pseudotime", y = NULL, title = "Signature score-based clusters") +
  theme_classic(base_size = 9) +
  theme(
    plot.title         = element_text(size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 9),
    legend.position    = "none"
  )
save_plot(p_clusters, "fig1H_clusters_vs_pseudotime.pdf", width = 10, height = 4)


# ==============================================================================
# Panel I - Genes upregulated in the Senescent cluster
# ==============================================================================

# --- Differential expression: Senescent vs all other cells ---
Idents(multiome_pt) <- "senescence_cluster"
senescent_markers <- FindMarkers(multiome_pt, ident.1 = "Senescent", only.pos = TRUE,
                                 min.pct = 0.1, logfc.threshold = 0.25)
senescent_markers$gene <- rownames(senescent_markers)
write.csv(senescent_markers, file.path(outdir, "suppl_senescent_markers.csv"), row.names = FALSE)

sig_senescent_genes <- senescent_markers$gene[
  senescent_markers$p_val_adj < 0.05 & senescent_markers$avg_log2FC > 2
]

# --- Overlap with the reference signatures (saved as a supplementary table) ---
# NB: 'X0.006149198' looks like an auto-generated column name from a headerless
# numeric column on import - confirm this is the intended filter.
universal_aging <- read.csv(path_universal_aging)
universal_aging <- universal_aging$Leng9[universal_aging$X0.006149198 > 0.01]

other_signatures <- unique(c(universal_aging, tabula_muris, sen_mayo, genes_up))
sig_senescent_overlap <- intersect(sig_senescent_genes, other_signatures)
write.csv(data.frame(gene = sig_senescent_overlap),
          file.path(outdir, "suppl_senescent_signature_overlap.csv"), row.names = FALSE)

# --- Dot plots by functional category ---
# The first factor level plots at the bottom, so "Senescent" is the bottom row.
multiome_pt$senescence_cluster <- factor(multiome_pt$senescence_cluster,
                                         levels = c("Senescent", "Acute injury", "Normal"))

sen_dotplot <- function(features, color, file, width) {
  missing_sig <- setdiff(features, sig_senescent_genes)
  if (length(missing_sig) > 0) {
    message(file, ": not in the Senescent marker set (padj < 0.05, log2FC > 2): ",
            paste(missing_sig, collapse = ", "))
  }
  p <- DotPlot(multiome_pt, group.by = "senescence_cluster", features = features,
               cols = c("grey90", color), scale.min = 0) +
    RotatedAxis() +
    guides(color = "none") +  # keep the dot-size legend, drop the colour bar
    theme(axis.text.x = element_text(face = "italic"),
          legend.title = element_blank())
  save_plot(p, file, width = width, height = 5)
}

sen_dotplot(
  c("TRAF1", "BIRC3", "NFKBID", "TNFAIP3", "KLF6", "CREB5", "HDAC9",          # NF-kB signalling
    "GFRA1", "DCLK1", "RAPGEF5", "SIPA1L2", "DUSP10", "ACSM1", "NCEH1",       # growth signalling
    "PLA2G4C", "OGFRL1"),
  "navy", "fig1I_dotplot_signalling.pdf", width = 16)

sen_dotplot(
  c("CXCL1", "CCL2", "IL18", "C3", "ICAM1", "GDF15", "KITLG", "TNFAIP2",      # SASP
    "TNFAIP8", "FAS", "HAVCR1"),
  "red4", "fig1I_dotplot_sasp.pdf", width = 12.5)

sen_dotplot(
  c("ITGA3", "ITGB8", "CLDN1", "LAMC2", "LAMA3", "COL12A1", "HS3ST1", "NRCAM"),  # ECM / adhesion
  "purple4", "fig1I_dotplot_ecm.pdf", width = 11)

sen_dotplot(
  c("ANXA2", "MAP1B", "TPM4", "FHL1", "MICAL2",                                # cytoskeletal
    "MAMLD1", "CELF2", "XKR4", "SLC35F3", "TMEM178B", "KIAA0408"),             # uncharacterised
  "darkgreen", "fig1I_dotplot_cytoskeleton_other.pdf", width = 13)


# ==============================================================================
# Panel J - CDKN1A chromatin accessibility (ATAC) and expression
# ==============================================================================

# --- Transfer cluster labels onto the full multiome object by barcode ---
labels <- setNames(rep(NA_character_, ncol(multiome)), colnames(multiome))
shared <- intersect(colnames(multiome_pt), colnames(multiome))
labels[shared] <- as.character(multiome_pt$senescence_cluster[shared])
multiome$senescence_cluster <- labels

multiome_labeled <- subset(multiome, cells = colnames(multiome)[!is.na(multiome$senescence_cluster)])

# --- Coverage plot (first level plots as the top track: Normal on top) ---
DefaultAssay(multiome_labeled) <- "ATAC"
multiome_labeled$senescence_cluster <- factor(multiome_labeled$senescence_cluster,
                                              levels = c("Normal", "Acute injury", "Senescent"))
Idents(multiome_labeled) <- "senescence_cluster"

p_atac <- CoveragePlot(
  multiome_labeled, region = "CDKN1A", group.by = "senescence_cluster",
  extend.upstream = 21000, extend.downstream = 1000,
  annotation = TRUE, peaks = FALSE, links = FALSE
) & scale_fill_manual(values = cluster_colors)
save_plot(p_atac, "fig1J_cdkn1a_atac.pdf", width = 10, height = 4.5)

# --- Companion dot plot (first level plots at the bottom: Senescent on the bottom) ---
multiome_labeled$senescence_cluster <- factor(multiome_labeled$senescence_cluster,
                                              levels = c("Senescent", "Acute injury", "Normal"))
p_cdkn1a <- DotPlot(multiome_labeled, features = "CDKN1A", group.by = "senescence_cluster",
                    assay = "SCT", cols = c("grey90", "navy"), scale.min = 0) +
  guides(color = "none")
save_plot(p_cdkn1a, "fig1J_cdkn1a_dotplot.pdf", width = 10, height = 4.5)


