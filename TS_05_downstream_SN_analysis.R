# STEP 1: Pseudobulk ----
set.seed(123)

suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(purrr); library(readr)
  library(SummarizedExperiment)
})

dir.create("PB", showWarnings = FALSE)

allcombined <- readRDS(file = "Data/all_cell_scores_new3.rds") # for CPP scores
path_milo <- "Data/milo_objects/v2_Immune.Mic.P2RY12.rds"   # Milo neighborhoods
path_npa  <- "Data/NPA/v2_Immune.Mic.P2RY12.rds"            # Differential abundance test results
path_seurat  <- "Data/mathys/Immune_cells.rds"              # Seurat

milo_obj <- readRDS(path_milo)
npa      <- readRDS(path_npa)
seu      <- readRDS(path_seurat)

milo_obj <- inject_umap_into_milo(milo_obj, seu) # function defined in TS_04

sn_dir <- "Data/superneighborhood_rds"

tab_both           <- readRDS(file.path(sn_dir, "tab_both.rds")) # mapping between neighborhoods and SNs
super_centers_umap <- readRDS(file.path(sn_dir, "super_centers_umap.rds"))
super_centers_pca  <- readRDS(file.path(sn_dir, "super_centers_pca.rds"))
cent_pca           <- readRDS(file.path(sn_dir, "cent_pca.rds"))
cent_umap          <- readRDS(file.path(sn_dir, "cent_umap.rds"))
da2                <- readRDS(file.path(sn_dir, "da_table_filtered.rds"))

# quick sanity checks
stopifnot(all(c("Nhood", "super_id", "set") %in% colnames(tab_both)))
stopifnot(all(c("UMAP_1", "UMAP_2", "super_id") %in% colnames(super_centers_umap)))
stopifnot("PCA_1" %in% colnames(super_centers_pca))

message("Loaded superneighborhood data from ", sn_dir)

# parameters
min_cells_per_donor_in_SN   <- 30
min_cells_per_donor_in_rest <- 30
min_donors_required         <- 10

# helpers
.map_nhood_to_cols <- function(nhood_ids, nh_mat) {
  cols  <- colnames(nh_mat)
  guess <- suppressWarnings(match(nhood_ids, as.integer(cols)))
  if (any(is.na(guess))) guess <- as.integer(nhood_ids)
  guess
}
cells_for_nhoods <- function(nhood_ids, nh_mat) {
  if (length(nhood_ids) == 0) return(character(0))
  col_idx <- .map_nhood_to_cols(nhood_ids, nh_mat)
  col_idx <- col_idx[is.finite(col_idx) & col_idx >= 1 & col_idx <= ncol(nh_mat)]
  if (!length(col_idx)) return(character(0))
  hit <- Matrix::rowSums(nh_mat[, col_idx, drop = FALSE] != 0) > 0
  rownames(nh_mat)[as.vector(hit)]
}

nh_mat <- milo_obj@nhoods # cell X neighborhood membership matrix
sn2nh  <- split(tab_both$Nhood, tab_both$super_id) # list of neighborhoods in each SN
sn2cells <- purrr::map(sn2nh, cells_for_nhoods, nh_mat = nh_mat) # list of cells in each SN

proj_by_cell <- as.character(SummarizedExperiment::colData(milo_obj)$projid) 
names(proj_by_cell) <- colnames(milo_obj)
cnts <- SummarizedExperiment::assay(milo_obj, "counts")

sum_cols <- function(cell_vec) Matrix::rowSums(cnts[, cell_vec, drop = FALSE])

sn_ids <- sort(names(sn2cells))
qc_rows <- list()

for (sn in sn_ids) {
  message("PB: ", sn)
  cells_in_sn <- sn2cells[[sn]]
  cells_rest  <- setdiff(colnames(milo_obj), cells_in_sn)
  donor_in_sn   <- split(cells_in_sn,  proj_by_cell[cells_in_sn])
  donor_in_rest <- split(cells_rest,   proj_by_cell[cells_rest])
  donors <- intersect(names(donor_in_sn), names(donor_in_rest))
  
  donors_keep <- donors[
    vapply(donor_in_sn[donors],   length, integer(1)) >= min_cells_per_donor_in_SN &
      vapply(donor_in_rest[donors], length, integer(1)) >= min_cells_per_donor_in_rest
  ]
  
  qc_rows[[sn]] <- tibble::tibble(
    super_id = sn,
    donor = donors,
    n_cells_sn   = vapply(donor_in_sn[donors],   length, integer(1)),
    n_cells_rest = vapply(donor_in_rest[donors], length, integer(1)),
    keep = donor %in% donors_keep
  )
  
  if (length(donors_keep) < min_donors_required) {
    message("  skipped (not enough paired donors).")
    next
  }
  
  sn_mat   <- do.call(cbind, lapply(donors_keep, \(d) sum_cols(donor_in_sn[[d]])))
  rest_mat <- do.call(cbind, lapply(donors_keep, \(d) sum_cols(donor_in_rest[[d]])))
  colnames(sn_mat)   <- paste0(sn, "__", donors_keep, "__SN")
  colnames(rest_mat) <- paste0(sn, "__", donors_keep, "__REST")
  
  pb <- list(
    counts  = cbind(sn_mat, rest_mat),
    samples = tibble::tibble(
      sample   = c(colnames(sn_mat), colnames(rest_mat)),
      donor    = rep(donors_keep, 2),
      group    = rep(c("SN","REST"), each = length(donors_keep)),
      super_id = sn
    )
  )
  saveRDS(pb, file = file.path("PB", paste0("PB_", sn, ".rds")))
}

qc <- dplyr::bind_rows(qc_rows)
readr::write_csv(qc, "PB/pseudobulk_qc_by_donor.csv")

# Additional QC on pseudobulks ----

suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr); library(ggplot2)
  library(purrr); library(readr); library(forcats); library(patchwork)
})

out_dir <- "PB"
dir.create(out_dir, showWarnings = FALSE)

# Highlighitng SNs of interest (based on inspection)
sn_set <- character(0) # not used

# Compute per-cell QC metrics
stopifnot(exists("milo_obj"))
cnts <- SummarizedExperiment::assay(milo_obj, "counts")
genes <- rownames(cnts)
cells <- colnames(cnts)

nCount   <- Matrix::colSums(cnts)
nFeature <- Matrix::colSums(cnts > 0)
is_mt    <- grepl("^MT-", genes, ignore.case = TRUE)
mt_counts <- if (any(is_mt)) Matrix::colSums(cnts[is_mt, , drop = FALSE]) else rep(0, length(cells))
percent_mt <- ifelse(nCount > 0, 100 * (mt_counts / nCount), NA_real_)

cell_qc <- tibble(
  cell         = cells,
  nCount_RNA   = as.numeric(nCount),
  nFeature_RNA = as.numeric(nFeature),
  percent_mt   = as.numeric(percent_mt)
)

# Map cells to super-neighborhoods
stopifnot(exists("tab_both"), "Nhood" %in% names(tab_both), "super_id" %in% names(tab_both))
nh_mat <- milo_obj@nhoods

.map_nhood_to_cols <- function(nhood_ids, nh_mat) {
  cols  <- colnames(nh_mat)
  guess <- suppressWarnings(match(nhood_ids, as.integer(cols)))
  if (any(is.na(guess))) guess <- as.integer(nhood_ids)
  guess
}
cells_for_nhoods <- function(nhood_ids, nh_mat) {
  if (length(nhood_ids) == 0) return(character(0))
  col_idx <- .map_nhood_to_cols(nhood_ids, nh_mat)
  col_idx <- col_idx[is.finite(col_idx) & col_idx >= 1 & col_idx <= ncol(nh_mat)]
  if (!length(col_idx)) return(character(0))
  hit <- Matrix::rowSums(nh_mat[, col_idx, drop = FALSE] != 0) > 0
  rownames(nh_mat)[as.vector(hit)]
}

sn2nh <- split(tab_both$Nhood, tab_both$super_id)
sn2cells <- purrr::map(sn2nh, cells_for_nhoods, nh_mat = nh_mat)

sn_long <- purrr::imap_dfr(sn2cells, ~tibble(cell = .x, super_id = .y))
sn_levels <- sort(unique(sn_long$super_id))
sn_long$super_id <- factor(sn_long$super_id, levels = sn_levels)

# Merge QC metrics with SN memberships
qc_by_sn <- sn_long %>%
  inner_join(cell_qc, by = "cell") %>%
  filter(is.finite(nCount_RNA), is.finite(nFeature_RNA)) %>%
  mutate(
    highlight = ifelse(super_id %in% sn_set, "highlight", "other"),
    highlight = factor(highlight, levels = c("other","highlight"))
  )

# plot
theme_qc <- theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

make_violin_highlight <- function(df, yvar, ylab, title = NULL, log10_y = FALSE, use_si = FALSE) {
  p <- ggplot(df, aes(x = super_id, y = .data[[yvar]], fill = highlight)) +
    geom_violin(trim = TRUE, color = "grey30", width = 0.8) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, color = "black") +
    scale_fill_manual(values = c("other" = "grey85", "highlight" = "#E64B35")) +
    labs(title = title, x = "Super-neighborhood", y = ylab) +
    theme_qc
  if (log10_y) {
    p <- p + scale_y_log10()
  } else if (use_si) {
    p <- p + scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")))
  }
  p
}

# build panels
p_nCount <- make_violin_highlight(
  qc_by_sn, "nCount_RNA", "nCount_RNA (UMIs)",
  title = "nCount_RNA by Super-neighborhood",
  log10_y = FALSE, use_si = TRUE
)

p_nFeature <- make_violin_highlight(
  qc_by_sn, "nFeature_RNA", "nFeature_RNA (# detected genes)",
  title = "nFeature_RNA by Super-neighborhood",
  log10_y = FALSE, use_si = FALSE
)

p_percent <- make_violin_highlight(
  qc_by_sn, "percent_mt", "percent.mt (%)",
  title = "percent.mt by Super-neighborhood",
  log10_y = FALSE, use_si = FALSE
)

# save
ggsave(file.path(out_dir, "qc_violin_highlight_nCount_by_SN.pdf"),   p_nCount,   width = 11, height = 5.5)
ggsave(file.path(out_dir, "qc_violin_highlight_nFeature_by_SN.pdf"), p_nFeature, width = 11, height = 5.5)
ggsave(file.path(out_dir, "qc_violin_highlight_percentMT_by_SN.pdf"), p_percent, width = 11, height = 5.5)

panel <- p_nCount / p_nFeature / p_percent + plot_layout(heights = c(1,1,1))
ggsave(file.path(out_dir, "qc_violin_highlight_panel_nCount_nFeature_percentMT_by_SN.pdf"),
       panel, width = 12, height = 12)

cat("\n SNs plotted:", paste(sn_set, collapse = ", "), "\nSaved in", out_dir, "\n")


# Build gene x SN matrix mat and gene z score matrix mat_scaled ----
set.seed(123)

suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(purrr)
  library(edgeR);  library(matrixStats)
  library(tidyr)
  library(ComplexHeatmap); library(circlize); library(grid)
})

pb_dir  <- "PB"  # directory containing PB_*.rds files created above
stopifnot(dir.exists(pb_dir))

# collect all pseudobulk files (all SNs)
pb_files <- list.files(pb_dir, pattern = "^PB_.*\\.rds$", full.names = TRUE)
stopifnot(length(pb_files) > 0)

# maintain stable SN order based on filenames
get_sn   <- function(f) sub("^PB_(.+)\\.rds$", "\\1", basename(f))
sn_order <- get_sn(pb_files)

# per-SN average logCPM profile
sn_profile_logcpm <- function(pb_rds) {
  pb <- readRDS(pb_rds)
  stopifnot(is.list(pb), !is.null(pb$counts), !is.null(pb$samples))
  sn <- unique(pb$samples$super_id); stopifnot(length(sn) == 1)
  sn_cols <- pb$samples$sample[pb$samples$group == "SN"]
  y <- edgeR::DGEList(pb$counts[, sn_cols, drop = FALSE])
  y <- edgeR::calcNormFactors(y, method = "TMM")
  lcp <- edgeR::cpm(y, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
  tibble::tibble(gene = rownames(lcp),
                 value = rowMeans(lcp),
                 super_id = sn)
}

profiles <- purrr::map(pb_files, sn_profile_logcpm)

# wide matrix across ALL SNs 
profiles_df <- dplyr::bind_rows(profiles) %>%
  dplyr::group_by(gene, super_id) %>%
  dplyr::summarise(value = mean(value), .groups = "drop")

profiles_df$super_id <- factor(profiles_df$super_id, levels = sn_order)

wide <- tidyr::pivot_wider(
  profiles_df,
  names_from  = super_id,
  values_from = value
)

# keep genes present in all SNs
wide <- wide[stats::complete.cases(wide), , drop = FALSE]

mat <- as.matrix(wide[, -1, drop = FALSE])
rownames(mat) <- wide$gene

# Z-score per gene (row-wise)
mat_scaled <- t(scale(t(mat)))
mat_scaled[!is.finite(mat_scaled)] <- 0

# STEP 2: Picking SNs for downstream analysis -----

set.seed(123)

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

stopifnot(exists("mat_scaled"))

## 1) Correlation matrix
S_cor <- cor(mat_scaled, method = "pearson", use = "pairwise.complete.obs")

## 2) Rename superneighborhoods
sn_names_original <- colnames(S_cor)

rename_sn <- function(x) {
  x <- gsub("^AD_SN_", "A", x)
  x <- gsub("^HC_SN_", "N", x)
  x
}

S_cor_disp <- S_cor
rownames(S_cor_disp) <- rename_sn(rownames(S_cor_disp))
colnames(S_cor_disp) <- rename_sn(colnames(S_cor_disp))

## 3) Robust color scale centered at 0
offdiag <- S_cor_disp[upper.tri(S_cor_disp) | lower.tri(S_cor_disp)]
rng <- stats::quantile(offdiag, c(0.02, 0.98), na.rm = TRUE)
col_fun_cor <- circlize::colorRamp2(
  c(rng[1], 0, rng[2]),
  c("#2c7fb8", "#f7f7f7", "#d7301f")
)

## 4) Text color helper
.pick_text_col <- function(val) {
  if (!is.finite(val)) return(NA_character_)
  if (abs(val) >= 0.6) "white" else "black"
}

## 5) Correlation heatmap with labeled cells
ht_cor <- Heatmap(
  S_cor_disp,
  name = "r",
  col = col_fun_cor,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_names_rot = 45,
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- S_cor_disp[i, j]
    if (is.finite(val)) {
      grid.text(sprintf("%.2f", val), x = x, y = y,
                gp = gpar(fontsize = 8, col = .pick_text_col(val)))
    }
  }
)

ht_cor_drawn <- draw(ht_cor)

# take the column order as indices and map back to original names
ord_idx <- ComplexHeatmap::column_order(ht_cor_drawn)
sn_order_corr <- sn_names_original[ord_idx]   


## 6) save
out_dir<-"Figures/Superneighborhoods/picking_SNs/"
pdf(file.path(out_dir,"cor_heatmap_with_values.pdf"), width = 8, height = 7); draw(ht_cor); dev.off()
png(file.path(out_dir,"cor_heatmap_with_values.png"), width = 2000, height = 1800, res = 300); draw(ht_cor); dev.off()

# PLOT: Boxplot of CPP scores within SNs ----
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(purrr)
})

# allcombined: rows=cells, columns: cpp scores
# sn2cells: named list super_id -> character vector of cells
# sn_order: character vector of super_id in heatmap order

make_cpp_boxplot <- function(allcombined, sn2cells, sn_order, cpp_col = "AD_status") {
  stopifnot(is.data.frame(allcombined), cpp_col %in% colnames(allcombined))
  stopifnot(is.list(sn2cells), length(sn2cells) > 0)
  stopifnot(is.character(sn_order), length(sn_order) > 0)
  
  # per-cell CPP vector
  cpp_vec <- allcombined[[cpp_col]]
  names(cpp_vec) <- rownames(allcombined)
  
  # align barcodes 
  align_cpp <- function(v, cells) {
    out <- v[cells]
    if (mean(is.na(out)) > 0.2) {
      strip <- function(x) sub("-\\d+-\\d+$", "", x)
      out <- v[match(strip(cells), strip(names(v)))]
      names(out) <- cells
    }
    out
  }
  
  # long per-cell table
  df <- purrr::imap_dfr(sn2cells, function(cells, sid) {
    if (!length(cells)) return(NULL)
    tibble(super_id = sid, cell = cells, cpp_val = align_cpp(cpp_vec, cells))
  }) %>% filter(is.finite(cpp_val))
  
  keep <- sn_order[sn_order %in% df$super_id]
  df$super_id <- factor(df$super_id, levels = rev(keep))
  
  ggplot(df, aes(x = super_id, y = cpp_val)) +
    
    geom_jitter(aes(color = cpp_val), width = 0.15, alpha = 0.1, size = 0.1, shape = 16) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "grey25",alpha=0) +
    scale_color_gradient2(low = "#0059FF", mid = "#E8E1DA", high = "#FF1818",
                          midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                          name = "Cell CPP") +
    geom_hline(yintercept = 0,linewidth = 0.5, color = "grey40") +
    ylim(-2,2)+
    coord_flip() +
    theme_classic() +
    theme(legend.position = "top") +
    labs(title = paste0("Per-cell CPP by super-neighborhood (", cpp_col, ")"),
         x = "Super-neighborhood (heatmap top→bottom)", y = paste0("CPP: ", cpp_col))
}

p <- make_cpp_boxplot(allcombined, sn2cells, sn_order = sn_order_corr, cpp_col = "AD_status")
ggsave("Figures/Superneighborhoods/CPP_boxplots_per_SN_AD_status.pdf", p, width = 2.4, height = 5.5)


# PLOT: Boxplot of cell donor origin composition of SNs
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

# load previous outputs
out_dir <- "Data/superneighborhood_rds"
sn_cell_with_donor <- readRDS(file.path(out_dir, "sn_cell_with_donor.rds"))

stopifnot(exists("combined_metadata"),
          all(c("projid","AD_status") %in% colnames(combined_metadata)))

# join projid + AD_status
sn_cell_with_donor <- sn_cell_with_donor %>%
  mutate(projid = as.integer(projid))

combined_md_min <- combined_metadata %>%
  transmute(projid = as.integer(projid),
            AD_status = as.numeric(AD_status))

sn_cells_ad <- sn_cell_with_donor %>%
  left_join(combined_md_min, by = "projid")

# Count cells per SN × AD_status
sn_counts <- sn_cells_ad %>%
  filter(!is.na(AD_status)) %>%
  group_by(super_id, AD_status) %>%
  summarise(n_cells = n(), .groups = "drop")

# Determine SN order
sn_counts <- sn_counts %>%
    mutate(super_id = factor(super_id, levels = rev(sn_order_corr)))  

# plot
pal_status <- c("0" = "#0D47A1", "1" = "#C21807")

p_counts <- ggplot(sn_counts, aes(x = super_id, y = n_cells, fill = factor(AD_status))) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = pal_status,
                    labels = c("Control (0)", "AD (1)"),
                    name = "AD status") +
  coord_flip() +
  scale_y_reverse() +  
  labs(x = "Super-neighborhood", y = "Number of cells",
       title = "Total cells per super-neighborhood (bars extend left)\ncolored by donor AD status") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))

print(p_counts)

# save
ggsave("Figures/Superneighborhoods/sn_cells_by_ADstatus.pdf", p_counts, width = 1.6, height = 5.5)

# PLOT: Jaccard Index heatmap of cell overlap among SNs ----
suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(ComplexHeatmap); library(circlize)
})

# Inputs
in_dir <- "Data/superneighborhood_rds"
memb <- readRDS(file.path(in_dir, "sn_cell_with_donor.rds")) %>%
  dplyr::select(super_id, cell) %>% distinct()
stopifnot(all(c("super_id","cell") %in% names(memb)))
stopifnot(exists("sn_order_corr"))

# Incidence matrix
sn_levels   <- sort(unique(memb$super_id))
cell_levels <- sort(unique(memb$cell))
M <- sparseMatrix(
  i = match(memb$cell, cell_levels),
  j = match(memb$super_id, sn_levels),
  x = 1,
  dims = c(length(cell_levels), length(sn_levels)),
  dimnames = list(cell_levels, sn_levels)
)

# Jaccard index
S <- as.matrix(crossprod(M))
n <- Matrix::colSums(M)
U <- outer(n, n, "+") - S
U[U == 0] <- NA_real_
J <- S / U
diag(J) <- 1

# sort by sn_order_corr
keep <- intersect(sn_order_corr, colnames(J))
J <- J[keep, keep, drop = FALSE]

# Heatmap 
col_fun <- colorRamp2(c(0, 1), c("white", "black"))
cell_size <- unit(0.4, "cm")

ht <- Heatmap(
  J,
  name = "Jaccard",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1)),
  width  = cell_size * ncol(J),
  height = cell_size * nrow(J)
)

# Save PDF
res_dir <- "Results/Superneighborhoods"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(res_dir, "SN_jaccard_heatmap.pdf"),
    width = convertWidth(cell_size * ncol(J), "in", valueOnly = TRUE) + 1.5,
    height = convertHeight(cell_size * nrow(J), "in", valueOnly = TRUE) + 1.5)
draw(ht)
dev.off()


# PLOT: Jaccard Index heatmap of donor overlap among SNs ----
suppressPackageStartupMessages({
  library(dplyr); library(Matrix); library(ComplexHeatmap); library(circlize); library(grid)
})

# input
in_dir  <- "Data/superneighborhood_rds"
res_dir <- "Figures/Superneighborhoods"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

memb <- readRDS(file.path(in_dir, "sn_cell_with_donor.rds")) %>%
  dplyr::select(super_id, projid) %>%
  dplyr::filter(!is.na(super_id), !is.na(projid)) %>%
  dplyr::distinct()

stopifnot(all(c("super_id", "projid") %in% names(memb)))

# Incidence matrix: projid × super_id (1 if SN contains ≥1 cell from that donor)
sn_levels   <- sort(unique(memb$super_id))
don_levels  <- sort(unique(memb$projid))

P <- sparseMatrix(
  i = match(memb$projid,  don_levels),
  j = match(memb$super_id, sn_levels),
  x = 1,
  dims = c(length(don_levels), length(sn_levels)),
  dimnames = list(don_levels, sn_levels)
)

# Jaccard on donor sets per SN
S <- as.matrix(crossprod(P))      # |A ∩ B| (shared donors)
n <- Matrix::colSums(P)           # |A| per SN (unique donors)
U <- outer(n, n, "+") - S         # |A ∪ B| = |A| + |B| - |A∩B|
U[U == 0] <- NA_real_             
J <- S / U
diag(J) <- 1

# order
keep <- intersect(sn_order_corr, colnames(J))
J <- J[keep, keep, drop = FALSE]


# heatmap
col_fun   <- colorRamp2(c(0, 1), c("white", "black"))
cell_size <- unit(0.40, "cm")

ht <- Heatmap(
  J,
  name = "Jaccard",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1)),
  width  = cell_size * ncol(J),
  height = cell_size * nrow(J),
  use_raster = FALSE
)

# save
outfile_pdf <- file.path(res_dir, "SN_donor_jaccard_heatmap.pdf")
pdf(
  outfile_pdf,
  width  = convertWidth(cell_size * ncol(J), "in", valueOnly = TRUE) + 1.5,
  height = convertHeight(cell_size * nrow(J), "in", valueOnly = TRUE) + 1.5
)
draw(ht)
dev.off()

# PLOT: UMAP of selected SNs ----
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(Matrix)
  library(dplyr); library(purrr); library(tidyr)
  library(ggplot2); library(glue)
  library(svglite)      
})

stopifnot(
  exists("milo_obj"),
  "UMAP" %in% SingleCellExperiment::reducedDimNames(milo_obj),
  !is.null(milo_obj@nhoods),
  exists("tab_both"),
  all(c("super_id","Nhood") %in% colnames(tab_both))
)

sn_select <- c("AD_SN_5","AD_SN_2","AD_SN_3","HC_SN_8","HC_SN_7","HC_SN_9")
sn_colors <- c(
  "AD_SN_5" = "#C59207",
  "AD_SN_2" = "#F88E2B",
  "AD_SN_3" = "#9B2226",
  "HC_SN_8" = "#004466",
  "HC_SN_7" = "#57B2C4",
  "HC_SN_9" = "#28A474"
)

pt_bg   <- 0.8
pt_sel  <- 0.8
alpha_bg  <- 0.5
alpha_sel <- 0.8

out_dir <- "Figures/Superneighborhoods/Selected_UMAP"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## HELPERS
# Clean labels
rename_sn <- function(x) {
  x <- gsub("^AD_SN_", "A", x)
  x <- gsub("^HC_SN_", "N", x)
  x
}

# Map Nhood IDs to columns of nh_mat 
.map_nhood_to_cols <- function(nhood_ids, nh_mat) {
  cols  <- colnames(nh_mat)
  guess <- suppressWarnings(match(nhood_ids, as.integer(cols)))
  if (any(is.na(guess))) guess <- as.integer(nhood_ids)
  guess
}

# Union of cells for a set of neighborhoods
cells_for_nhoods <- function(nhood_ids, nh_mat) {
  if (length(nhood_ids) == 0) return(character(0))
  col_idx <- .map_nhood_to_cols(nhood_ids, nh_mat)
  col_idx <- col_idx[is.finite(col_idx) & col_idx >= 1 & col_idx <= ncol(nh_mat)]
  if (!length(col_idx)) return(character(0))
  hit <- Matrix::rowSums(nh_mat[, col_idx, drop = FALSE] != 0) > 0
  rownames(nh_mat)[as.vector(hit)]
}

## DATA
# UMAP for all cells (background)
emb_all <- as.data.frame(SingleCellExperiment::reducedDim(milo_obj, "UMAP")) %>%
  setNames(c("UMAP_1","UMAP_2")) %>%
  mutate(cell = colnames(milo_obj)) %>%
  filter(is.finite(UMAP_1), is.finite(UMAP_2))

# Neighborhood incidence matrix
nh_mat <- milo_obj@nhoods
stopifnot(!is.null(nh_mat), nrow(nh_mat) == ncol(milo_obj))

# Ensure selected SNs exist
all_sids <- unique(tab_both$super_id)
missing  <- setdiff(sn_select, all_sids)
if (length(missing)) stop("These superneighborhoods are not present in tab_both: ",
                          paste(missing, collapse = ", "))

# assign colors for SNs
if (!all(sn_select %in% names(sn_colors))) {
  left <- setdiff(sn_select, names(sn_colors))
  auto_cols <- scales::hue_pal()(length(left))
  names(auto_cols) <- left
  sn_colors <- c(sn_colors, auto_cols)
}

# Build cell membership for each selected SN
sn2cells <- map(sn_select, function(sid) {
  nh_in_sn <- tab_both %>% filter(super_id == sid) %>% pull(Nhood) %>% unique()
  cells_for_nhoods(nh_in_sn, nh_mat)
})
names(sn2cells) <- sn_select

# Resolve overlaps by priority (earlier SNs in sn_select take precedence)
cell_to_sn <- list()
for (sid in sn_select) {
  for (c in sn2cells[[sid]]) {
    if (is.null(cell_to_sn[[c]])) cell_to_sn[[c]] <- sid
  }
}
cell_to_sn <- unlist(cell_to_sn, use.names = TRUE)

sel_df <- emb_all %>%
  filter(cell %in% names(cell_to_sn)) %>%
  mutate(super_id = unname(cell_to_sn[cell]),
         super_lbl = rename_sn(super_id))

# Factor order in legend follows sn_select
sel_df$super_id  <- factor(sel_df$super_id, levels = sn_select)
sel_df$super_lbl <- factor(sel_df$super_lbl, levels = rename_sn(sn_select))

# Legend colors keyed to *labels* shown (clean names)
lab_colors <- sn_colors[sn_select]
names(lab_colors) <- rename_sn(names(lab_colors))

# Background without selected cells (prevents halos)
bg_df <- emb_all %>% filter(!cell %in% sel_df$cell)

## PLOT 
p <- ggplot() +
  # background (grey)
  geom_point(
    data = bg_df,
    aes(UMAP_1, UMAP_2),
    shape = 16, size = pt_bg, alpha = alpha_bg, color = "grey50"
  ) +
  # selected SN cells (colored)
  geom_point(
    data = sel_df,
    aes(UMAP_1, UMAP_2, color = super_lbl),
    shape = 16, size = pt_sel, alpha = alpha_sel
  ) +
  scale_color_manual(values = lab_colors, name = "Superneighborhood") +
  coord_fixed() +
  theme_void(base_size = 11) +
  theme(legend.position = "top")

tag <- paste(rename_sn(sn_select), collapse = "_")
pdf_file <- file.path(out_dir, glue("UMAP_SN_{tag}.pdf"))
svg_file <- file.path(out_dir, glue("UMAP_SN_{tag}.svg"))

# save
ggsave(pdf_file, p, width = 8.5, height = 8.5, device = cairo_pdf)

message("Saved outputs: ", pdf_file, " and ", svg_file)
p

# STEP 3: Microglial States ----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(fs)
  library(seriation)
})

stopifnot(exists("mat"), is.matrix(mat), ncol(mat) >= 2)
stopifnot(exists("sn_order_corr"))

# --- Display label mapping 
rename_sn <- function(x) {
  x <- gsub("^AD_SN_", "A", x)
  x <- gsub("^HC_SN_", "N", x)
  x
}

# ---- setup
set.seed(123)
path_markers <- "Data/microglia_states_Sun_2023.csv"
out_dir      <- "Results/MicrogliaStates"
out_pdf_blk  <- file.path(out_dir, "heatmap_markers_by_SN_diagStates.pdf")

low_col  <- "#2c7fb8"
mid_col  <- "#f7f7f7"
high_col <- "#d7301f"

# microglia state levels + palette
state_levels <- c("MG0","MG1","MG2","MG3","MG4","MG5","MG6",
                  "MG7","MG8","MG10","MG11","MG12")
state_palette <- c(
  "MG0"  = "#E22027","MG1"  = "#387FBB","MG2"  = "#4FAE49","MG3"  = "#99509F",
  "MG4"  = "#F57F20","MG5"  = "#F6EC3B","MG6"  = "#A85726","MG7"  = "#ED84B5",
  "MG8"  = "#989898","MG10" = "#6AC3A5","MG11" = "#F68E66","MG12" = "#8E9FCA"
)

# labelled marker genes
genes_to_label  <- c(
  "P2RY12","CX3CR1","FRMD4A","INO80D","TMEM163","CPEB4",
  "FTL","FTH1","MYO1E","PTPRG","CD163","F13A1","HSP90AA1",
  "HSPH1","NAMPT","SLC2A3","SPON1","LRRK2","CCL3","IL1B",
  "IFI44L1","MX1","EZH2","BRIP1"
)
label_with_state <- FALSE

# ---- Column order & display labels (identical to correlation heatmap) 
sn_order     <- sn_order_corr
sn_order_use <- sn_order[sn_order %in% colnames(mat)]
if (!length(sn_order_use)) stop("No overlap between sn_order_corr and colnames(mat).")
col_labs <- rename_sn(sn_order_use)

# ---- Load markers and build marker matrix from mat
state_markers <- readr::read_csv(path_markers, show_col_types = FALSE) %>%
  dplyr::select(microgliaState, gene, avg_log2FC) %>%
  dplyr::filter(microgliaState %in% state_levels) %>%
  dplyr::distinct(microgliaState, gene, .keep_all = TRUE) %>%
  dplyr::mutate(microgliaState = factor(microgliaState, levels = state_levels)) %>%
  dplyr::arrange(microgliaState, dplyr::desc(avg_log2FC), gene)

sel <- state_markers$gene %in% rownames(mat)
plot_rows <- state_markers[sel, , drop = FALSE]
if (!nrow(plot_rows)) stop("None of the marker genes are in 'mat'.")

# allow duplicated genes across states
row_keys <- make.unique(paste0(plot_rows$gene, " | ", as.character(plot_rows$microgliaState)))

expr_mark <- mat[plot_rows$gene, sn_order_use, drop = FALSE]
rownames(expr_mark) <- row_keys

# Row z-score within marker set
Z <- t(scale(t(expr_mark)))
Z[!is.finite(Z)] <- 0

# state factor per row (before reordering)
row_states <- factor(plot_rows$microgliaState, levels = state_levels)
names(row_states) <- rownames(Z)

# Diagonalization helpers

order.tsp <- function(mat, rows = TRUE, method = "euclidean") {
  if (rows) {
    tsp <- seriation::seriate(dist(mat, method = method))
  } else {
    tsp <- seriation::seriate(dist(t(mat), method = method))
  }
  ord <- seriation::get_order(tsp, 1)
  if (rows) {
    return(mat[ord, , drop = FALSE])
  } else {
    return(mat[, ord, drop = FALSE])
  }
}

diag.mat3 <- function(mat, ratio = 0.5, cutoff = 0.25, rows = TRUE) {
  prev.cutoff <- attr(mat, "cutoff")
  if (rows) {
    ord <- order(
      rowSums(mat > cutoff) > ratio * ncol(mat),
      apply(mat, 1, which.max),
      decreasing = TRUE
    )
    mat <- mat[ord, , drop = FALSE]
    cto <- apply(mat, 1, which.max)
    idx <- rowSums(mat > cutoff) > ratio * ncol(mat)
    cto[idx] <- 0
  } else {
    ord <- order(
      colSums(mat > cutoff) > ratio * nrow(mat),
      apply(mat, 2, which.max),
      decreasing = TRUE
    )
    mat <- mat[, ord, drop = FALSE]
    cto <- apply(mat, 2, which.max)
    idx <- colSums(mat > cutoff) > ratio * nrow(mat)
    cto[idx] <- 0
  }
  attr(mat, "cutoff") <- prev.cutoff
  if (is.null(attr(mat, "cutoff"))) {
    if (rows) {
      attr(mat, "cutoff") <- list(row = cto, col = NULL)
    } else {
      attr(mat, "cutoff") <- list(row = NULL, col = cto)
    }
  } else {
    if (rows) {
      attr(mat, "cutoff")$row <- cto
    } else {
      attr(mat, "cutoff")$col <- cto
    }
  }
  return(mat)
}

htSortMatrix <- function(M, method = "euclidean", ratio = 0.5, cutoff = 0.25, sort = c(1, 2)) {
  if (is.logical(sort)) {
    if (sort) {
      sort <- c(1, 2)
    } else {
      sort <- c()
    }
  }
  rows    <- 1 %in% sort
  columns <- 2 %in% sort
  if (rows) {
    M <- order.tsp(M, rows = TRUE, method = method)
  }
  if (columns) {
    M <- order.tsp(M, rows = FALSE, method = method)
  }
  if (rows) {
    M <- diag.mat3(M, rows = TRUE, ratio = ratio, cutoff = cutoff)
  }
  if (columns) {
    M <- diag.mat3(M, rows = FALSE, ratio = ratio, cutoff = cutoff)
  }
  return(M)
}

# Diagonalize states (keeping individual state blocks here and just sorting those)

# states that have markers in Z
states_present <- intersect(state_levels, unique(as.character(row_states)))

# state-level summary: mean z across genes in each state, per SN
if (length(states_present) < 2L) {
  warning("Fewer than 2 states present; skipping diagonalization of state blocks.")
  state_order <- states_present
} else {
  S <- vapply(
    states_present,
    function(st) colMeans(Z[row_states == st, , drop = FALSE]),
    numeric(ncol(Z))
  )
  S <- t(S)  # rows = states, cols = SNs
  rownames(S) <- states_present
  
  # diagonalize rows (states)
  S_sorted <- htSortMatrix(S, method = "euclidean", ratio = 0.5, cutoff = 0.25, sort = 1)
  state_order <- rownames(S_sorted)
}

# reorder rows of Z by new state block order, keeping genes within a state together
ord_rows <- order(factor(row_states, levels = state_order))
Z_blk    <- Z[ord_rows, , drop = FALSE]
row_states_blk <- row_states[ord_rows]

# flip the row order after diagonalization
Z_blk          <- Z_blk[nrow(Z_blk):1, , drop = FALSE]
row_states_blk <- rev(row_states_blk)

# refactor with new (flipped) block order for splitting/legend
row_states_blk <- factor(
  as.character(row_states_blk),
  levels = rev(state_order)  # flip the state block order as well
)

# Labels in the new row order (after flip)
annot_idx_blk <- which(gsub(" \\| .*$", "", rownames(Z_blk)) %in% genes_to_label)
labels_to_show_blk <- if (label_with_state) {
  paste0(
    gsub(" \\| .*$", "", rownames(Z_blk)[annot_idx_blk]),
    " (", as.character(row_states_blk[annot_idx_blk]), ")"
  )
} else {
  gsub(" \\| .*$", "", rownames(Z_blk)[annot_idx_blk])
}


# Heatmap

fs::dir_create(out_dir, recurse = TRUE)

col_fun <- circlize::colorRamp2(
  c(-1.5, 0, 1.5),
  c(low_col, mid_col, high_col)
)

ra_state_blk <- rowAnnotation(
  State = row_states_blk,
  col = list(State = state_palette),
  gp  = gpar(col = NA),
  annotation_legend_param = list(title = "Microglia state")
)

ra_marks_blk <- rowAnnotation(
  Gene = anno_mark(
    at         = annot_idx_blk,
    labels     = labels_to_show_blk,
    which      = "row",
    side       = "right",
    labels_gp  = gpar(fontsize = 6),
    link_width = unit(6, "mm"),
    link_gp    = gpar(lwd = 0.6, col = "#444444"),
    extend     = unit(c(0, 2), "mm")
  )
)

ht_core_blk <- Heatmap(
  Z_blk,
  name              = "row-z",
  col               = col_fun,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  column_order      = sn_order_use,   # matches correlation heatmap
  column_labels     = col_labs,
  show_row_names    = FALSE,
  show_column_names = TRUE,
  column_names_gp   = gpar(fontsize = 8),
  row_split         = row_states_blk,
  heatmap_legend_param = list(
    title  = "z",
    at     = c(-1.5, -1, 0, 1, 1.5),
    labels = c("-1.5","-1","0","1","1.5")
  ),
  use_raster = FALSE
)

ht_blk <- ra_state_blk + ht_core_blk + ra_marks_blk

draw(ht_blk)

# PDF
grDevices::pdf(out_pdf_blk, width = 8, height = 6)
draw(ht_blk)
grDevices::dev.off()

cat(sprintf(
  "Block-diagonal states heatmap: %d SNs, %d marker rows, %d states\nOrder of states: %s\nOutput: %s\n",
  ncol(Z_blk), nrow(Z_blk), length(state_order),
  paste(state_order, collapse = ", "),
  out_pdf_blk
))
################################################################################
################################################################################
# STEP 4: DEGs ----
################################################################################
################################################################################


# setup
source("000 utils.r")
prepare_metadata <- function(phenotypes = c("AD_status","amyloid","tangles","plaq_n","cogn_global_lv")) {
  meta <- load_meta_data()
  revised <- load_celltype_revised_phenotype(phenotypes)
  revised_unified <- do.call(
    cbind,
    lapply(names(revised), function(nm) {
      df <- revised[[nm]]
      colnames(df) <- paste(nm, colnames(df), sep = ".")
      df
    })
  )
  
  meta$AD_status <- ifelse(meta$niareagansc %in% c(4,3),0,1) ## Binarize AD status
  meta <- meta %>% mutate(projid = as.integer(projid))
  rownames(meta) <- meta$projid
  combined_metadata <- meta %>%
    select(projid, AD_status, amyloid, tangles, plaq_n, cogn_global_lv, pmi, msex, age_death) %>%
    inner_join(
      revised_unified %>% tibble::rownames_to_column("projid") %>% mutate(projid = as.integer(projid)),
      by = "projid"
    )
  rownames(combined_metadata) <- combined_metadata$projid
  return(combined_metadata)
}

phenotypes <- c("AD_status","amyloid","tangles","plaq_n","cogn_global_lv")
combined_metadata <- prepare_metadata(phenotypes)

# DEG analysis
suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr); library(purrr); library(tibble); library(jsonlite)
})

dir.create("Results/SN_DEG/qc",  recursive = TRUE, showWarnings = FALSE)
dir.create("Results/SN_DEG/logs", recursive = TRUE, showWarnings = FALSE)

pairs <- list(
  c("AD_SN_5","HC_SN_8"),
  c("AD_SN_2","HC_SN_7"),
  c("AD_SN_3","HC_SN_9")
)

load_pb <- function(sid) {
  pb <- readRDS(file.path("PB", paste0("PB_", sid, ".rds")))
  stopifnot("counts" %in% names(pb), "samples" %in% names(pb))
  idx <- pb$samples$group == "SN"
  counts <- pb$counts[, pb$samples$sample[idx], drop = FALSE]
  meta   <- pb$samples[idx, c("sample","donor","super_id"), drop = FALSE]
  list(counts = counts, meta = meta)
}

run_one <- function(A, B, combined_metadata) {
  # ---- helpers ----
  dedup_rows <- function(mat) {
    # collapse duplicated gene IDs by summing counts
    if (anyDuplicated(rownames(mat))) {
      mat <- rowsum(mat, group = rownames(mat))
    }
    mat
  }
  sanitize_meta <- function(cm) {
    cm <- as.data.frame(cm)
    # Ensure exactly one projid column (integer)
    if (!"projid" %in% colnames(cm)) {
      cm <- tibble::rownames_to_column(cm, var = "projid")
    }
    cm$projid <- as.integer(cm$projid)
    
    # Keep only needed covariates to avoid accidental duplicates
    need <- c("projid", "pmi", "msex", "age_death")
    ok   <- intersect(need, colnames(cm))
    cm <- cm[, ok, drop = FALSE]
    
    # Name-repair (just in case upstream supplied dupes)
    tibble::as_tibble(cm, .name_repair = "unique")
  }
  
  load_pb <- function(sid) {
    pb <- readRDS(file.path("PB", paste0("PB_", sid, ".rds")))
    stopifnot("counts" %in% names(pb), "samples" %in% names(pb))
    idx <- pb$samples$group == "SN"
    counts <- pb$counts[, pb$samples$sample[idx], drop = FALSE]
    counts <- dedup_rows(counts)
    meta   <- pb$samples[idx, c("sample","donor","super_id"), drop = FALSE]
    list(counts = counts, meta = meta)
  }
  
  # ---- load & harmonize ----
  pbA <- load_pb(A); pbB <- load_pb(B)
  
  genes <- union(rownames(pbA$counts), rownames(pbB$counts))
  cA <- matrix(0, nrow = length(genes), ncol = ncol(pbA$counts),
               dimnames = list(genes, colnames(pbA$counts)))
  cB <- matrix(0, nrow = length(genes), ncol = ncol(pbB$counts),
               dimnames = list(genes, colnames(pbB$counts)))
  cA[rownames(pbA$counts), colnames(pbA$counts)] <- pbA$counts
  cB[rownames(pbB$counts), colnames(pbB$counts)] <- pbB$counts
  counts <- cbind(cA, cB)
  
  cm <- sanitize_meta(combined_metadata)
  
  # Build colData WITHOUT creating duplicate names
  colData <- dplyr::bind_rows(
    dplyr::mutate(pbA$meta, group = A),
    dplyr::mutate(pbB$meta, group = B)
  ) |>
    dplyr::transmute(
      sample,
      projid = as.integer(donor),
      group
    ) |>
    dplyr::left_join(cm, by = "projid") |>
    dplyr::mutate(group = factor(group, levels = c(B, A))) |>
    tibble::column_to_rownames("sample")
  
  # Synchronize columns with counts
  counts <- counts[, rownames(colData), drop = FALSE]
  
  # Drop samples with missing covariates
  keep <- stats::complete.cases(colData[, c("pmi","msex","age_death")])
  counts <- counts[, keep, drop = FALSE]
  colData <- colData[keep, , drop = FALSE]
  
  nA <- sum(colData$group == A); nB <- sum(colData$group == B)
  if (nA < 3 || nB < 3) {
    qc <- list(skipped = TRUE, reason = "too_few_donors", nA = nA, nB = nB,
               comparison = paste0(A,"_vs_",B))
    return(list(res = NULL, qc = qc))
  }
  
  # ---- DESeq2 ----
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData   = colData,
    design    = ~ group + pmi + msex + age_death
  )
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("group", A, B)) |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      comparison = paste0(A,"_vs_",B),
      n_donors_groupA = nA,
      n_donors_groupB = nB
    )
  
  qc <- list(
    comparison = paste0(A,"_vs_",B),
    n_genes = nrow(dds),
    n_samples = ncol(dds),
    n_donors_groupA = nA,
    n_donors_groupB = nB,
    design = "~ group + pmi + msex + age_death"
  )
  list(res = res, qc = qc)
}

all_res <- list()
for (pr in pairs) {
  A <- pr[1]; B <- pr[2]
  out_base <- file.path("Results/SN_DEG", paste0("DESeq2_", A, "_vs_", B))
  out_qc   <- file.path("Results/SN_DEG/qc", paste0("DESeq2_", A, "_vs_", B, "_runmeta.json"))
  
  ans <- run_one(A, B, combined_metadata)
  if (is.null(ans$res)) {
    writeLines(jsonlite::toJSON(ans$qc, auto_unbox = TRUE, pretty = TRUE), out_qc)
    next
  }
  saveRDS(ans$res, paste0(out_base, ".rds"))
  write_csv(ans$res, paste0(out_base, ".csv"))
  writeLines(jsonlite::toJSON(ans$qc, auto_unbox = TRUE, pretty = TRUE), out_qc)
  all_res[[paste0(A,"_vs_",B)]] <- ans$res
}

if (length(all_res)) {
  combo <- bind_rows(all_res)
  saveRDS(combo, "Results/SN_DEG/DESeq2_ALL_SNpairs.rds")
  write_csv(combo, "Results/SN_DEG/DESeq2_ALL_SNpairs.csv")
}

# Volcano Plots
suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(dplyr); library(purrr)
  library(ggrastr)  
})

harmonize_deseq <- function(df) {
  tiny <- .Machine$double.xmin
  out <- df %>%
    transmute(
      gene   = if ("gene" %in% names(df)) gene else rownames(df),
      log2FC = as.numeric(if ("log2FC" %in% names(df)) log2FC else log2FoldChange),
      FDR    = as.numeric(if ("FDR"    %in% names(df)) FDR    else padj)
    )
  out$mlog10FDR <- -log10(pmax(out$FDR, tiny))
  out
}

compute_global_limits <- function(all_res_list, fdr_thr = 0.01, x_pad = 0.25) {
  dd <- bind_rows(lapply(all_res_list, harmonize_deseq)) %>%
    filter(is.finite(log2FC), is.finite(mlog10FDR), is.finite(FDR))
  if (!nrow(dd)) return(list(x = c(-2,2), y = c(0,5)))
  xmax <- max(abs(dd$log2FC), na.rm = TRUE); xmax <- ceiling((xmax + x_pad) * 2) / 2
  ymax <- max(dd$mlog10FDR, na.rm = TRUE);  ymax <- ceiling(ymax)
  list(x = c(-xmax, xmax), y = c(0, ymax))
}

make_volcano <- function(df,
                         title,
                         colors         = list(pos = "#D55E00", neg = "#0072B2", ns = "grey70"),
                         point_size     = 1.6,
                         point_alpha    = 0.85,
                         xlim_fixed,
                         ylim_fixed,
                         fdr_thr        = 0.01,
                         label_per_side = 20,
                         fdr_fill2      = 0.20,
                         seed           = 123) {
  
  stopifnot(is.data.frame(df))
  
  tiny <- .Machine$double.xmin
  d <- transform(
    data.frame(df, check.names = FALSE),
    gene   = if ("gene" %in% names(df)) df[["gene"]] else
      if (!is.null(rownames(df))) rownames(df) else NA_character_,
    log2FC = as.numeric(if ("log2FC" %in% names(df)) df[["log2FC"]] else df[["log2FoldChange"]]),
    FDR    = as.numeric(if ("FDR"    %in% names(df)) df[["FDR"]]    else df[["padj"]])
  )
  
  # keep only finite FDR and log2FC; compute y as -log10(FDR)
  d <- d[is.finite(d$FDR) & is.finite(d$log2FC), , drop = FALSE]
  d$mlog10FDR <- -log10(pmax(d$FDR, tiny))
  
  # significance direction (FDR only)
  d$sig_dir <- "ns"
  d$sig_dir[d$FDR < fdr_thr & d$log2FC > 0] <- "pos_sig"
  d$sig_dir[d$FDR < fdr_thr & d$log2FC < 0] <- "neg_sig"
  pal <- c(pos_sig = colors$pos, neg_sig = colors$neg, ns = colors$ns)
  
  # ---- plot ----
  p <- ggplot2::ggplot(d, ggplot2::aes(x = log2FC, y = mlog10FDR)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    ggplot2::geom_vline(xintercept = 0,               linetype = "solid",  linewidth = 0.5, color = "grey40") +
    ggrastr::geom_point_rast(ggplot2::aes(color = sig_dir),
                             alpha = point_alpha, size = point_size,
                             stroke = 0, dpi = 600, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::coord_cartesian(xlim = xlim_fixed, ylim = ylim_fixed, clip = "off") +
    ggplot2::labs(title = title, x = expression(log[2](FC)), y = expression(-log[10](FDR))) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0))
  
  pick_side <- function(side_df, n) {
    ord <- function(z) z[order(z$FDR, -abs(z$log2FC)), , drop = FALSE]
    t1 <- subset(side_df, FDR < fdr_thr)
    t2 <- subset(side_df, FDR < fdr_fill2)
    t3 <- side_df
    out <- ord(t1)
    if (nrow(out) < n) out <- unique(rbind(out, ord(t2)))
    if (nrow(out) < n) out <- unique(rbind(out, ord(t3)))
    head(out, n)
  }
  
  pool  <- d[!is.na(d$gene) & d$gene != "", , drop = FALSE]
  up_df <- pick_side(subset(pool, log2FC > 0),  label_per_side)
  dn_df <- pick_side(subset(pool, log2FC < 0),  label_per_side)
  lab_df <- rbind(up_df, dn_df)
  lab_df$clean_gene <- gsub("\\.", " ", lab_df$gene)
  
  if (nrow(lab_df)) {
    set.seed(seed)
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      ggplot2::aes(label = clean_gene, color = sig_dir),
      size = 3, max.overlaps = 100, box.padding = 0.35, point.padding = 0.25,
      force = 0.9, min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6
    )
  }
  
  n_up   <- sum(d$sig_dir == "pos_sig", na.rm = TRUE)
  n_down <- sum(d$sig_dir == "neg_sig", na.rm = TRUE)
  message(sprintf("[volcano] %s — FDR<%.2f: up=%d, down=%d; labeled up=%d, down=%d",
                  title, fdr_thr, n_up, n_down, nrow(up_df), nrow(dn_df)))
  
  return(p)
}


fdr_thr <- 0.01
lims <- compute_global_limits(all_res, fdr_thr = fdr_thr)

plot_dir <- "Results/SN_DEG/volcano_deseq"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

p1 <- make_volcano(
  all_res[["AD_SN_5_vs_HC_SN_8"]],
  title = "AD_SN_5 vs HC_SN_8",
  colors = list(pos = "#C2922D", neg = "#004668", ns = "grey70"),
  xlim_fixed = lims$x, ylim_fixed = lims$y,
  fdr_thr = fdr_thr, label_per_side = 20, fdr_fill2 = 0.20
)
ggsave(file.path(plot_dir, "volcano_AD_SN_5_vs_HC_SN_8.pdf"),
       p1, width = 6, height = 5, device = cairo_pdf)

p2 <- make_volcano(
  all_res[["AD_SN_2_vs_HC_SN_7"]],
  title = "AD_SN_2 vs HC_SN_7",
  colors = list(pos = "#F78E2B", neg = "#57B1C3", ns = "grey70"),
  xlim_fixed = lims$x, ylim_fixed = lims$y,
  fdr_thr = fdr_thr, label_per_side = 20, fdr_fill2 = 0.20
)
ggsave(file.path(plot_dir, "volcano_AD_SN_2_vs_HC_SN_7.pdf"),
       p2, width = 6, height = 5, device = cairo_pdf)

p3 <- make_volcano(
  all_res[["AD_SN_3_vs_HC_SN_9"]],
  title = "AD_SN_3 vs HC_SN_9",
  colors = list(pos = "#9C2328", neg = "#29A576", ns = "grey70"),
  xlim_fixed = lims$x, ylim_fixed = lims$y,
  fdr_thr = fdr_thr, label_per_side = 20, fdr_fill2 = 0.20
)
ggsave(file.path(plot_dir, "volcano_AD_SN_3_vs_HC_SN_9.pdf"),
       p3, width = 6, height = 5, device = cairo_pdf)

pdf(file.path(plot_dir, "volcanoes.gridArrange.pdf"), height=5, width=18)
grid.arrange(p1, p2,p3, nrow=1)
dev.off()
### counting DEGs

count_de_genes <- function(df, fdr_thr = 0.01) {
  # Handle column name variations (same as make_volcano)
  log2FC <- as.numeric(if ("log2FC" %in% names(df)) df[["log2FC"]] else df[["log2FoldChange"]])
  FDR    <- as.numeric(if ("FDR"    %in% names(df)) df[["FDR"]]    else df[["padj"]])
  
  # Keep only finite values
  keep <- is.finite(FDR) & is.finite(log2FC)
  FDR <- FDR[keep]
  log2FC <- log2FC[keep]
  
  # Count by direction
  n_up   <- sum(FDR < fdr_thr & log2FC > 0, na.rm = TRUE)
  n_down <- sum(FDR < fdr_thr & log2FC < 0, na.rm = TRUE)
  
  data.frame(
    up = n_up,
    down = n_down,
    total_de = n_up + n_down
  )
}

count_de_genes(all_res[["AD_SN_5_vs_HC_SN_8"]])
count_de_genes(all_res[["AD_SN_2_vs_HC_SN_7"]])
count_de_genes(all_res[["AD_SN_3_vs_HC_SN_9"]])

###### comparing DEG sets ######
library(eulerr)
library(ggvenn)

get_de_genes <- function(df, fdr_thr = 0.01, direction = "both") {
  genes <- if ("gene" %in% names(df)) df[["gene"]] else
    if (!is.null(rownames(df))) rownames(df) else NULL
  log2FC <- as.numeric(if ("log2FC" %in% names(df)) df[["log2FC"]] else df[["log2FoldChange"]])
  FDR    <- as.numeric(if ("FDR"    %in% names(df)) df[["FDR"]]    else df[["padj"]])
  keep <- is.finite(FDR) & is.finite(log2FC) & !is.na(genes)
  genes <- genes[keep]
  FDR <- FDR[keep]
  log2FC <- log2FC[keep]
  sig <- FDR < fdr_thr
  if (direction == "up") {
    sig <- sig & log2FC > 0
  } else if (direction == "down") {
    sig <- sig & log2FC < 0
  }
  
  return(genes[sig])
}

deg_sets <- list(
  "AD_SN_5_vs_HC_SN_8" = get_de_genes(all_res[["AD_SN_5_vs_HC_SN_8"]]),
  "AD_SN_2_vs_HC_SN_7" = get_de_genes(all_res[["AD_SN_2_vs_HC_SN_7"]]),
  "AD_SN_3_vs_HC_SN_9" = get_de_genes(all_res[["AD_SN_3_vs_HC_SN_9"]])
)

deg_up <- list(
  "Homeostatic-like" = get_de_genes(all_res[["AD_SN_5_vs_HC_SN_8"]], direction = "up"),
  "Inflammatory-like" = get_de_genes(all_res[["AD_SN_2_vs_HC_SN_7"]], direction = "up"),
  "Glycolytic-like" = get_de_genes(all_res[["AD_SN_3_vs_HC_SN_9"]], direction = "up")
)

deg_down <- list(
  "Homeostatic-like" = get_de_genes(all_res[["AD_SN_5_vs_HC_SN_8"]], direction = "down"),
  "Inflammatory-like" = get_de_genes(all_res[["AD_SN_2_vs_HC_SN_7"]], direction = "down"),
  "Glycolytic-like" = get_de_genes(all_res[["AD_SN_3_vs_HC_SN_9"]], direction = "down")
)

# Print sizes to verify
cat("Up in AD-like SNs:\n")
sapply(deg_up, length)

cat("\nUp in Non-AD-like SNs:\n")
sapply(deg_down, length)

# Create plots with clearer styling
p1 <- ggvenn(deg_up, 
             fill_color = c("#C59207", "#F88E2B", "#9B2226"),
             stroke_size = 0.75, 
             set_name_size = 5,
             text_size = 4) + 
  ggtitle("Up-regulated in AD-like SNs") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

p2 <- ggvenn(deg_down, 
             fill_color = c("#004466", "#57B2C4", "#28A474"),
             stroke_size = 0.75, 
             set_name_size = 5,
             text_size = 4) + 
  ggtitle("Up-regulated in Non-AD-like SNs") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

pdf("pdf/sn.degOverlaps.pdf", height=6, width=8)
grid.arrange(p1, p2, nrow = 1)
dev.off()

################################################################################
################################################################################
# STEP 5: FGSEA on Hallmark pathways ----
################################################################################
################################################################################
suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(readr); library(tibble)
  library(msigdbr); library(fgsea); library(ggplot2)
})

set.seed(123)

# setup 
minSize <- 5L; maxSize <- 1000L
top_n_per_panel <- 10
low_col  <- "grey80"; high_col <- "black"

out_dir   <- "Results/SN_DEG/fgsea"
plot_file <- file.path(out_dir, "fgsea_faceted_by_contrast.pdf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ranks
make_ranks <- function(deseq, gene_col = "gene",
                       lfc_col = "log2FoldChange", p_col = "pvalue") {
  df <- as.data.frame(deseq)
  stopifnot(all(c(gene_col, lfc_col, p_col) %in% names(df)))
  genes <- as.character(df[[gene_col]])
  lfc   <- suppressWarnings(as.numeric(df[[lfc_col]]))
  pval  <- suppressWarnings(as.numeric(df[[p_col]]))
  tiny  <- .Machine$double.xmin
  pval[!is.finite(pval) | is.na(pval) | pval <= 0] <- 1
  base  <- sign(lfc) * (-log10(pmax(pval, tiny)))
  eps   <- 1e-9 * rank(abs(lfc), ties.method = "first")
  score <- base + sign(lfc) * eps
  keep  <- is.finite(score) & nzchar(genes)
  ranks <- score[keep]; names(ranks) <- genes[keep]
  sort(ranks, decreasing = TRUE)
}

# fgsea runner
run_fgsea_on_res <- function(res, label, gsets,
                             minSize = 5, maxSize = 1000) {
  ranks   <- make_ranks(res)
  all_ids <- unique(unlist(gsets, use.names = FALSE))
  ranks   <- ranks[names(ranks) %in% all_ids]
  if (!length(ranks)) return(tibble::tibble())
  fg <- fgsea::fgseaMultilevel(
    pathways    = gsets,
    stats       = ranks,
    minSize     = minSize,
    maxSize     = maxSize
  )
  dplyr::as_tibble(fg) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES))) |>
    dplyr::mutate(contrast = label,
                  direction = ifelse(NES >= 0, "up", "down"))
}

# make gene sets
as_pathway_list <- function(df, id_col = "gene_symbol") {
  stopifnot(all(c("gs_id", id_col) %in% names(df)))
  tmp <- unique(df[, c("gs_id", id_col)])
  sp  <- split(tmp[[id_col]], tmp[["gs_id"]])
  lapply(sp, function(v) unique(v[!is.na(v) & nzchar(v)]))
}

species <- "Homo sapiens"

df_H <- msigdbr(species = species, collection = "H") |>
  transmute(collection = "H", gs_name, gene_symbol, gs_id = paste0("H:", gs_name))

df_GOBP <- msigdbr(species = species, collection = "C5", subcollection = "BP") |>
  transmute(collection = "GOBP", gs_name, gene_symbol, gs_id = paste0("GOBP:", gs_name))

gset_df <- bind_rows(df_H, df_GOBP)

gsets_combined <- as_pathway_list(gset_df)

gs_map <- gset_df |> distinct(gs_id, collection, gs_name)

# run fgsea
fgsea_list <- purrr::imap(all_res, function(df, contrast_name) {
  res <- run_fgsea_on_res(df, contrast_name, gsets_combined,
                          minSize = minSize, maxSize = maxSize)
  if (!nrow(res)) return(NULL)
  # 'pathway' returned by fgsea is our gs_id; join to recover collection/name
  res |>
    mutate(gs_id = .data$pathway) |>
    left_join(gs_map, by = "gs_id") |>
    # Use a readable label while keeping origin visible
    mutate(pathway_display = paste0("[", collection, "] ", gs_name))
})

fgsea_tbl <- dplyr::bind_rows(fgsea_list)
stopifnot(nrow(fgsea_tbl) > 0)

# Save tables (keep origin info)
readr::write_csv(
  fgsea_tbl |>
    dplyr::select(contrast, collection, gs_name, everything(), -leadingEdge, -pathway_display),
  file.path(out_dir, "fgsea_ALL_contrasts_H_plus_GOBP.csv")
)

readr::write_csv(
  fgsea_tbl |>
    filter(padj < 0.05) |>
    arrange(contrast, padj, desc(abs(NES))) |>
    dplyr::select(contrast, collection, gs_name, NES, padj, size, direction),
  file.path(out_dir, "fgsea_SIG_padj_lt_0.05_H_plus_GOBP.csv")
)

# Main Plot

sn_colors <- c(
  "AD_SN_5" = "#C59207",
  "AD_SN_2" = "#F88E2B",
  "AD_SN_3" = "#9B2226",
  "HC_SN_8" = "#004466",
  "HC_SN_7" = "#57B2C4",
  "HC_SN_9" = "#28A474"
)

parse_contrast <- function(x) {
  ab <- strsplit(x, "_vs_", fixed = TRUE)[[1]]
  if (length(ab) != 2) ab <- c(NA_character_, NA_character_)
  tibble(contrast = x, A = ab[1], B = ab[2])
}

# Top-N per contrast (combined collections)
plot_tbl <- fgsea_tbl %>%
  group_by(contrast) %>%
  arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = top_n_per_panel) %>%
  ungroup()

contrast_order <- names(all_res)
plot_tbl <- plot_tbl %>%
  mutate(
    contrast = factor(contrast, levels = contrast_order),
    is_sig     = padj < 0.05,
    mlog10padj = -log10(padj)
  )

# A/B map + fill color by NES direction (same logic as before)
ab_map <- unique(as.character(plot_tbl$contrast)) %>%
  lapply(parse_contrast) %>% dplyr::bind_rows()

plot_tbl <- plot_tbl %>%
  left_join(ab_map, by = "contrast") %>%
  mutate(fill_col = ifelse(NES >= 0, sn_colors[A], sn_colors[B]))

# Per-contrast x-axis ordering (top → bottom by padj then |NES|)
plot_tbl <- plot_tbl %>%
  group_by(contrast) %>%
  arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
  mutate(
    x_key = paste0(pathway_display, "|||", contrast),
    x_key = factor(x_key, levels = rev(unique(x_key)))
  ) %>%
  ungroup()

# Layers to show non-sig outlined and sig filled by A/B color with alpha by FDR
plot_sig <- filter(plot_tbl, is_sig)
plot_all <- plot_tbl

max_abs_nes <- max(abs(plot_all$NES), na.rm = TRUE)
nes_lim <- c(-max_abs_nes, max_abs_nes)

strip_suffix <- function(x) sub("\\|\\|\\|.*$", "", x)

p <- ggplot(plot_all, aes(x = x_key, y = NES)) +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4, colour = "grey40") +
  geom_col(fill = "white", colour = "black", linewidth = 0.25) +
  geom_col(data = plot_sig, aes(fill = fill_col, alpha = mlog10padj), colour = NA) +
  coord_flip() +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
  scale_y_continuous(limits = nes_lim) +
  scale_fill_identity(guide = "none") +
  scale_alpha_continuous(range = c(0.25, 1), name = expression(-log[10](FDR))) +
  scale_x_discrete(labels = strip_suffix) +
  labs(x = "Pathway", y = "NES",
       title = "fgsea per contrast (Hallmark + GO:BP combined)\nnon-sig = white+outline; sig = SN color (A if NES>0, B if NES<0), alpha = −log10(FDR)") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    panel.spacing = unit(0.6, "lines")
  )

ggsave(plot_file, p, width = 15, height = 10, device = cairo_pdf)
message("Wrote plot: ", plot_file)


# PLOT: Distribution of cells across selected super-neighborhoods by AD status ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

# Inputs
out_dir <- "Data/superneighborhood_rds"
sn_cell_with_donor <- readRDS(file.path(out_dir, "sn_cell_with_donor.rds"))
stopifnot(exists("combined_metadata"),
          all(c("projid","AD_status") %in% colnames(combined_metadata)))

# Choose SNs and colors
sn_select <- c("AD_SN_5","AD_SN_2","AD_SN_3","HC_SN_8","HC_SN_7","HC_SN_9")
sn_colors <- c(
  "AD_SN_5" = "#C59207",
  "AD_SN_2" = "#F88E2B",
  "AD_SN_3" = "#9B2226",
  "HC_SN_8" = "#004466",
  "HC_SN_7" = "#57B2C4",
  "HC_SN_9" = "#28A474"
)

# Join donor AD status
df0 <- sn_cell_with_donor %>%
  mutate(projid = as.integer(projid)) %>%
  left_join(combined_metadata %>%
              transmute(projid = as.integer(projid),
                        AD_status = as.numeric(AD_status)),
            by = "projid") %>%
  filter(AD_status %in% c(0,1),
         super_id %in% sn_select) %>%
  dplyr::select(AD_status, super_id)

# Count and normalize within AD_status
counts <- df0 %>%
  count(AD_status, super_id, name = "n_cells") %>%
  complete(AD_status, super_id = sn_select, fill = list(n_cells = 0)) %>%
  group_by(AD_status) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(AD_status_f = ifelse(AD_status == 1, "AD donors", "Control donors"),
         super_id = factor(super_id, levels = sn_select))

# Plot two pies (AD vs Control)
p <- ggplot(counts, aes(x = "", y = prop, fill = super_id)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ AD_status_f, nrow = 1) +
  scale_fill_manual(values = sn_colors, name = "Super-neighborhood") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold"))

# Save 
res_dir <- "Figures/Superneighborhoods"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(res_dir, "sn_membership_pies_selected.pdf"),
       p, width = 9, height = 5)


################################################################################
# PLOT: UMAP grid of all SNs ----
################################################################################
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(glue)
  library(cowplot)
  library(tibble)
})

# config

cpp_col <- "AD_status"  

pt_bg     <- 0.3
pt_sel    <- 0.4
alpha_bg  <- 0.3
alpha_sel <- 0.9

out_dir <- "Figures/Superneighborhoods/SN_UMAP_CPP"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

max_cells_cloud <- 200000   # subsample background for speed

# checks

stopifnot(
  exists("milo_obj"),
  "UMAP" %in% SingleCellExperiment::reducedDimNames(milo_obj),
  !is.null(milo_obj@nhoods),
  exists("tab_both"),
  all(c("super_id","Nhood") %in% colnames(tab_both)),
  exists("allcombined"),
  exists("sn_order_corr"),
  is.character(sn_order_corr)
)

# helpers

# Clean labels like before (for titles)
rename_sn <- function(x) {
  x <- gsub("^AD_SN_", "A", x)
  x <- gsub("^HC_SN_", "N", x)
  x
}

# Map Nhood IDs to columns of nh_mat 
.map_nhood_to_cols <- function(nhood_ids, nh_mat) {
  cols  <- colnames(nh_mat)
  guess <- suppressWarnings(match(nhood_ids, as.integer(cols)))
  if (any(is.na(guess))) guess <- as.integer(nhood_ids)
  guess
}

# Union of cells for a set of neighborhoods
cells_for_nhoods <- function(nhood_ids, nh_mat) {
  if (length(nhood_ids) == 0) return(character(0))
  col_idx <- .map_nhood_to_cols(nhood_ids, nh_mat)
  col_idx <- col_idx[is.finite(col_idx) & col_idx >= 1 & col_idx <= ncol(nh_mat)]
  if (!length(col_idx)) return(character(0))
  hit <- Matrix::rowSums(nh_mat[, col_idx, drop = FALSE] != 0) > 0
  rownames(nh_mat)[as.vector(hit)]
}


# data

# UMAP for all cells
emb_all <- as.data.frame(SingleCellExperiment::reducedDim(milo_obj, "UMAP")) %>%
  setNames(c("UMAP_1","UMAP_2")) %>%
  mutate(cell = colnames(milo_obj)) %>%
  filter(is.finite(UMAP_1), is.finite(UMAP_2))

# Neighborhood incidence matrix
nh_mat <- milo_obj@nhoods
stopifnot(!is.null(nh_mat), nrow(nh_mat) == ncol(milo_obj))

# superneighborhoods that  exist
all_sids <- unique(tab_both$super_id)
sn_ids   <- intersect(sn_order_corr, all_sids)
if (!length(sn_ids)) {
  stop("No overlap between sn_order_corr and tab_both$super_id.")
}

ac <- as.data.frame(allcombined)
if (!"cell" %in% colnames(ac)) {
  ac <- tibble::rownames_to_column(ac, "cell")
}
stopifnot(
  "cell" %in% colnames(ac),
  cpp_col %in% colnames(ac)
)

cell_cpp <- emb_all %>%
  left_join(ac %>% dplyr::select(cell, !!cpp_col), by = "cell") %>%
  rename(cpp_val = !!cpp_col)

# Symmetric limits over all cells (for comparable colors across SNs)
rng <- range(cell_cpp$cpp_val, na.rm = TRUE)
cpp_max <- max(abs(rng))
if (!is.finite(cpp_max) || cpp_max == 0) {
  cpp_limits <- c(-1, 1)
} else {
  cpp_limits <- c(-cpp_max, cpp_max)
}

cpp_scale_common <- scale_color_gradient2(
  low = "#0059FF",
  mid = "#E8E1DA",
  high = "#FF1818",
  midpoint = 0,
  limits = cpp_limits,
  oob = scales::squish,
  name = cpp_col
)

emb_bg <- emb_all
if (nrow(emb_bg) > max_cells_cloud) {
  set.seed(123)
  emb_bg <- emb_bg[sample.int(nrow(emb_bg), max_cells_cloud), ]
}


make_sn_umap <- function(sid) {
  sid_chr <- as.character(sid)
  
  nh_in_sn <- tab_both %>%
    filter(super_id == sid_chr) %>%
    pull(Nhood) %>%
    unique()
  if (!length(nh_in_sn)) return(NULL)
  
  cells_in_sn <- cells_for_nhoods(nh_in_sn, nh_mat)
  if (!length(cells_in_sn)) return(NULL)
  
  cells_sn_cpp <- cell_cpp %>%
    filter(cell %in% cells_in_sn)
  
  n_cells <- nrow(cells_sn_cpp)
  cap_txt <- glue("{sid_chr} • cells: {n_cells}")
  
  ggplot() +
    # background (all cells, grey)
    geom_point(
      data = emb_bg,
      aes(UMAP_1, UMAP_2),
      shape = 16,
      size  = pt_bg,
      alpha = alpha_bg,
      color = "grey70"
    ) +
    # SN cells, colored by CPP
    geom_point(
      data = cells_sn_cpp,
      aes(UMAP_1, UMAP_2, color = cpp_val),
      shape = 16,
      size  = pt_sel,
      alpha = alpha_sel
    ) +
    cpp_scale_common +
    coord_fixed() +
    theme_void(base_size = 11) +
    theme(
      legend.position = "none",           # no legends
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.caption    = element_text(hjust = 0.5)
    ) +
    labs(
      title   = rename_sn(sid_chr),
      caption = cap_txt
    )
}

# build pdf

plots <- lapply(sn_ids, make_sn_umap)
names(plots) <- sn_ids
plots <- plots[!vapply(plots, is.null, logical(1))]

if (!length(plots)) {
  stop("No plots created – check that sn_ids have cells.")
}

pdf_file <- file.path(out_dir, "SN_UMAP_CPP_1col_sn_order_corr.pdf")

page_width  <- 7
page_height <- 1.8 * length(plots)   

pdf(pdf_file, width = page_width, height = page_height)
cowplot::plot_grid(
  plotlist = plots,
  ncol     = 1,
  align    = "v"
)
dev.off()

message("Saved 1-column CPP-colored SN UMAP PDF (no legends) to: ", pdf_file)

