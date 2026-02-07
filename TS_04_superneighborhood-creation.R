#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggrastr)
  library(ggnewscale)
  library(Matrix)
  library(igraph)
  library(readr)
})

# ======= Functions =======================================
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


# ======= Parameters ======================================
# Inputs
path_milo <- "Data/milo_objects/v2_Immune.Mic.P2RY12.rds"   # Milo (SCE) with PCA + UMAP
path_npa  <- "Data/NPA/v2_Immune.Mic.P2RY12.rds"            # list of DA frames
path_seurat  <- "Data/mathys/Immune_cells.rds"              # Seurat
test_id   <- "Test_1"                                       # which DA result to use
source("000 utils.r")
phenotypes <- c("AD_status","amyloid","tangles","plaq_n","cogn_global_lv")
combined_metadata <- prepare_metadata(phenotypes)
stopifnot(all(c("projid","AD_status") %in% colnames(combined_metadata)))


# Selection of top-N neighborhoods per tail after FDR
fdr_thresh     <- 0.05
n_top_each     <- 1000

# Graph & clustering (using PCA space of Milo)
k_centroids    <- 24
min_super_size <- 50
leiden_res     <- 0.0092  # resolution parameter for Leiden

# Plotting
max_cells_cloud <- 200000
seed_use        <- 123
pal_set         <- c("AD-like" = "#C21807", "healthy-like" = "#0D47A1")

dir.create("Figures/Superneighborhoods", recursive = TRUE, showWarnings = FALSE)
set.seed(seed_use)

# ======= Helpers ==========================================
# Getting neighborhood centroids in any reduction (mean of member cells)
nhood_centroids_df <- function(milo, reduction = "PCA", ndims = NULL) {
  stopifnot(reduction %in% SingleCellExperiment::reducedDimNames(milo))
  emb <- SingleCellExperiment::reducedDim(milo, reduction)
  if (!is.null(ndims)) emb <- emb[, ndims, drop = FALSE]
  nh  <- milo@nhoods                         # rows = cells, cols = neighborhoods
  sums  <- Matrix::t(nh) %*% emb             # (nhoods × dims)
  sizes <- Matrix::colSums(nh)               # length = nhoods
  cents <- sweep(as.matrix(sums), 1, sizes, "/")
  out <- as.data.frame(cents)
  colnames(out) <- paste0(reduction, "_", seq_len(ncol(out)))
  out$Nhood <- seq_len(nrow(out))
  out
}

# Build undirected KNN graph on rows of a matrix (Euclidean)
build_knn <- function(coords, k = 10) {
  n <- nrow(coords)
  if (n < 2) return(igraph::make_empty_graph(n, directed = FALSE))
  D <- as.matrix(dist(coords, method = "euclidean"))
  diag(D) <- Inf
  nbrs <- lapply(seq_len(n), function(i) {
    j <- order(D[i, ])[seq_len(min(k, n - 1))]
    cbind(from = i, to = j, w = 1/(1 + D[i, j]))
  })
  E <- do.call(rbind, nbrs)
  E2 <- as.data.frame(E)
  key <- paste(pmin(E2$from, E2$to), pmax(E2$from, E2$to), sep = "_")
  E2 <- E2[!duplicated(key), , drop = FALSE]
  igraph::graph_from_data_frame(E2, directed = FALSE,
                                vertices = data.frame(id = seq_len(n)))
}

# Leiden clustering of the subsets of neighborhoods using PCA centroids
cluster_subset_pca <- function(nhood_ids, pca_centers, k = 12,
                               min_super_size = 5, prefix = "SN",
                               leiden_res = 1.0) {
  
  if (!length(nhood_ids)) {
    return(dplyr::tibble(Nhood = integer(0), super_id = character(0)))
  }
  
  H <- pca_centers %>%
    dplyr::filter(Nhood %in% nhood_ids) %>%
    dplyr::arrange(Nhood)
  
  coords <- as.matrix(H %>% dplyr::select(dplyr::starts_with("PCA_")))
  g   <- build_knn(coords, k = k)
  
  com <- igraph::cluster_leiden(
    g,
    weights = igraph::E(g)$w,
    resolution_parameter = leiden_res
  )
  memb <- igraph::membership(com)
  
  out  <- tibble::tibble(
    Nhood = H$Nhood,
    super_id = paste0(prefix, "_", as.integer(memb))
  )
  
  keep <- names(which(table(out$super_id) >= min_super_size))
  out  <- out %>% dplyr::filter(super_id %in% keep)
  
  out
}

# Map a vector of neighborhood IDs to column indices of milo_obj@nhoods
.map_nhood_to_cols <- function(nhood_ids, nh_mat) {
  cols  <- colnames(nh_mat)
  guess <- suppressWarnings(match(nhood_ids, as.integer(cols)))
  if (any(is.na(guess))) guess <- as.integer(nhood_ids)
  guess
}

# Collect the set of cell barcodes belonging to a union of neighborhoods
cells_for_nhoods <- function(nhood_ids, nh_mat) {
  if (length(nhood_ids) == 0) return(character(0))
  col_idx <- .map_nhood_to_cols(nhood_ids, nh_mat)
  col_idx <- col_idx[is.finite(col_idx) & col_idx >= 1 & col_idx <= ncol(nh_mat)]
  if (!length(col_idx)) return(character(0))
  hit <- Matrix::rowSums(nh_mat[, col_idx, drop = FALSE] != 0) > 0
  rownames(nh_mat)[as.vector(hit)]
}

# copy UMAP from Seurat into Milo
inject_umap_into_milo <- function(milo, seu) {
  if ("UMAP" %in% SingleCellExperiment::reducedDimNames(milo)) return(milo)
  stopifnot("umap" %in% names(seu@reductions))
  um <- Seurat::Embeddings(seu, "umap")[, 1:2, drop = FALSE]
  colnames(um) <- c("UMAP_1", "UMAP_2")
  
  # align names
  cn_milo <- colnames(milo); cn_seur <- rownames(um)
  idx <- match(cn_milo, cn_seur)
  if (anyNA(idx)) {
    strip <- function(x) sub("-\\d+-\\d+$", "", x)
    rownames(um) <- strip(rownames(um))
    idx <- match(strip(cn_milo), rownames(um))
  }
  
  um_aligned <- matrix(NA_real_, nrow = length(cn_milo), ncol = 2,
                       dimnames = list(cn_milo, c("UMAP_1","UMAP_2")))
  ok <- !is.na(idx)
  um_aligned[ok, ] <- um[idx[ok], , drop = FALSE]
  SingleCellExperiment::reducedDim(milo, "UMAP") <- um_aligned
  milo
}


# ======= Data ========================================
allcombined <- readRDS(file = "Data/all_cell_scores_new3.rds") # for CPP scores
milo_obj <- readRDS(path_milo)
npa      <- readRDS(path_npa)
seu      <- readRDS(path_seurat)

milo_obj <- inject_umap_into_milo(milo_obj, seu)

#checks
stopifnot(all(c("PCA","UMAP") %in% reducedDimNames(milo_obj)))
stopifnot(test_id %in% names(npa))
da <- npa[[test_id]]
stopifnot(all(c("Nhood","logFC","FDR") %in% colnames(da)))

# use all available pcs
n_pcs_avail <- ncol(reducedDim(milo_obj, "PCA"))
ndims_use   <- seq_len(n_pcs_avail)

# ======= Neighborhood centroids (PCA for logic, UMAP for viz) =======
cent_pca <- nhood_centroids_df(milo_obj, reduction = "PCA",  ndims = ndims_use) %>%
  dplyr::mutate(Nhood = as.integer(Nhood))
cent_umap <- nhood_centroids_df(milo_obj, reduction = "UMAP", ndims = 1:2) %>%
  dplyr::rename(UMAP_1 = UMAP_1, UMAP_2 = UMAP_2) %>%
  dplyr::mutate(Nhood = as.integer(Nhood))

# Attach DA to UMAP centroids 
da2 <- da %>% dplyr::mutate(Nhood = as.integer(Nhood))
cent_all <- cent_umap %>%
  dplyr::left_join(da2 %>% dplyr::select(Nhood, logFC, FDR), by = "Nhood")

# ======= Select AD-like & healthy-like (top-N after FDR) ============
sig_da <- da2 %>% dplyr::filter(is.finite(FDR), FDR < fdr_thresh, is.finite(logFC))

pos_ids <- sig_da %>% dplyr::arrange(dplyr::desc(logFC)) %>%
  dplyr::slice_head(n = n_top_each) %>% dplyr::pull(Nhood)
neg_ids <- sig_da %>% dplyr::arrange(logFC) %>%
  dplyr::slice_head(n = n_top_each) %>% dplyr::pull(Nhood)

tab_AD <- cent_all %>% dplyr::filter(Nhood %in% pos_ids) %>% dplyr::mutate(set = "AD-like")
tab_HC <- cent_all %>% dplyr::filter(Nhood %in% neg_ids) %>% dplyr::mutate(set = "healthy-like")
if (!nrow(tab_AD) && !nrow(tab_HC)) stop("No neighborhoods selected for either set.")

# ======= CLUSTER IN PCA SPACE (super-neighborhoods) =================
# PCA centroids aligned to selected neighborhoods
pca_AD <- cent_pca %>% dplyr::filter(Nhood %in% tab_AD$Nhood)
pca_HC <- cent_pca %>% dplyr::filter(Nhood %in% tab_HC$Nhood)

mem_AD <- cluster_subset_pca(tab_AD$Nhood, cent_pca, k = k_centroids,
                             min_super_size = min_super_size, prefix = "AD_SN",
                             leiden_res = leiden_res)
mem_HC <- cluster_subset_pca(tab_HC$Nhood, cent_pca, k = k_centroids,
                             min_super_size = min_super_size, prefix = "HC_SN",
                             leiden_res = leiden_res)

tab_AD2 <- tab_AD %>% dplyr::inner_join(mem_AD, by = "Nhood")
tab_HC2 <- tab_HC %>% dplyr::inner_join(mem_HC, by = "Nhood")

tab_both <- dplyr::bind_rows(tab_AD2, tab_HC2)
if (!nrow(tab_both)) stop("After size filtering, no super-neighborhoods remain.")
tab_both$set <- factor(tab_both$set, levels = c("AD-like","healthy-like"))

# --- Super centroids in PCA (mean PCA-centroid per SN)

super_centers_pca <- tab_both %>%
  dplyr::select(super_id, Nhood) %>%
  dplyr::left_join(cent_pca, by = "Nhood") %>%
  dplyr::group_by(super_id) %>%
  dplyr::summarize(dplyr::across(starts_with("PCA_"), ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop")

# For plotting: mean UMAP of member neighborhoods 
super_centers_umap <- tab_both %>%
  dplyr::group_by(super_id, set) %>%
  dplyr::summarize(
    UMAP_1    = mean(UMAP_1, na.rm = TRUE),
    UMAP_2    = mean(UMAP_2, na.rm = TRUE),
    super_lfc = stats::median(logFC, na.rm = TRUE),
    .groups   = "drop"
  )
super_centers_umap$set <- factor(super_centers_umap$set, levels = c("AD-like","healthy-like"))

# ======= Summarize cell membership by donor projid ======
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tibble); library(purrr)
  library(SummarizedExperiment)
})

# Per-cell metadata. Row names = cell barcodes.
cell_md <- as.data.frame(SummarizedExperiment::colData(milo_obj)) %>%
  tibble::rownames_to_column("cell")

stopifnot("projid" %in% colnames(cell_md))  


# For each super_id, list member cells and summarize by projid
sn_ids <- unique(tab_both$super_id)

sn_cells_list <- setNames(vector("list", length(sn_ids)), sn_ids)

for (sid in sn_ids) {
  nh <- tab_both %>% filter(super_id == sid) %>% pull(Nhood)
  sn_cells_list[[sid]] <- cells_for_nhoods(nh, milo_obj@nhoods)
}

# Long table: cell ↔ super_id
sn_cell_tbl <- imap_dfr(sn_cells_list, ~tibble(super_id = .y, cell = .x))

# Join projid and per-cell attributes from colData
sn_cell_with_donor <- sn_cell_tbl %>%
  left_join(cell_md %>% select(cell, projid), by = "cell")

# Per-SN × projid counts and fractions
sn_donor_summary <- sn_cell_with_donor %>%
  filter(!is.na(projid)) %>%
  group_by(super_id, projid) %>%
  summarise(n_cells = n(), .groups = "drop_last") %>%
  mutate(total_in_SN = sum(n_cells),
         frac        = n_cells / total_in_SN) %>%
  ungroup() %>%
  arrange(super_id, desc(n_cells))



# ======= SAVE MAIN OUTPUTS =====================================
out_dir <- "Data/superneighborhood_rds"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(tab_both,            file = file.path(out_dir, "tab_both.rds"))
saveRDS(super_centers_umap,  file = file.path(out_dir, "super_centers_umap.rds"))
saveRDS(super_centers_pca,   file = file.path(out_dir, "super_centers_pca.rds"))
saveRDS(cent_pca,            file = file.path(out_dir, "cent_pca.rds"))
saveRDS(cent_umap,           file = file.path(out_dir, "cent_umap.rds"))
saveRDS(da2,                 file = file.path(out_dir, "da_table_filtered.rds"))
saveRDS(sn_cell_with_donor,  file= file.path(out_dir,"sn_cell_with_donor.rds"))
saveRDS(sn_donor_summary,    file= file.path(out_dir,"sn_donor_summary.rds"))

message("Saved super-neighborhood data to ", out_dir)

# ======= PLOT: UMAP with PCA-defined SNs for each individual SN =====
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment)
  library(Matrix); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggrepel); library(ggrastr); library(ggnewscale)
  library(glue); library(readr)
})
# --------------------------- required inputs ---------------------------
stopifnot(exists("milo_obj"),
          "UMAP" %in% SingleCellExperiment::reducedDimNames(milo_obj),
          !is.null(milo_obj@nhoods),
          exists("tab_both"),
          all(c("super_id","Nhood","UMAP_1","UMAP_2","logFC","set") %in% colnames(tab_both)),
          exists("da"),
          all(c("Nhood","logFC") %in% colnames(da)))

# If not present, build centroids and super_centers_umap
if (!exists("cent_umap")) {
  nhood_centroids_df <- function(milo, reduction = "UMAP", ndims = 1:2) {
    emb <- SingleCellExperiment::reducedDim(milo, reduction)[, ndims, drop = TRUE]
    nh  <- milo@nhoods
    sums  <- Matrix::t(nh) %*% as.matrix(emb)
    sizes <- Matrix::colSums(nh)
    cents <- sweep(as.matrix(sums), 1, sizes, "/")
    out <- as.data.frame(cents)
    colnames(out) <- paste0(reduction, "_", seq_len(ncol(out)))
    out$Nhood <- seq_len(nrow(out))
    out
  }
  cent_umap <- nhood_centroids_df(milo_obj, "UMAP", 1:2) %>% mutate(Nhood = as.integer(Nhood))
}

if (!exists("super_centers_umap")) {
  stopifnot(all(c("super_id","Nhood","UMAP_1","UMAP_2","logFC","set") %in% colnames(tab_both)))
  super_centers_umap <- tab_both %>%
    group_by(super_id, set) %>%
    summarise(
      UMAP_1 = mean(UMAP_1, na.rm = TRUE),
      UMAP_2 = mean(UMAP_2, na.rm = TRUE),
      super_lfc = stats::median(logFC, na.rm = TRUE),
      .groups = "drop"
    )
}

if (!exists("pal_set")) pal_set <- c("AD-like" = "#C21807", "healthy-like" = "#0D47A1")

# --------------------------- output dirs --------------------------------
out_dir <- "Figures/Superneighborhoods/SN_2D_layers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- knobs --------------------------------------
cells_bg_max         <- 200000   # cap background cells
cells_per_neigh_max  <- 250      # per-neighborhood cap for cell→neigh lines
cells_total_cap      <- 15000    # global cap per SN for cell→neigh lines
show_cell_lines      <- TRUE
show_cell_points     <- TRUE
show_neigh_lines     <- TRUE

point_sizes <- list(
  cell_bg   = 0.02,
  cell_sn   = 0.15,
  neigh_bg  = 0.5,
  neigh_sn  = 1.6,
  super     = 4.5
)
alphas <- list(
  cell_bg   = 0.10,
  cell_sn   = 0.80,
  cell_line = 0.06,
  neigh_bg  = 0.20,
  neigh_sn  = 0.90,
  neigh_line= 0.50,
  super     = 0.95
)
strokes <- list(
  neigh_sn  = 0.25,
  super     = 1.1
)

bwr_fill <- scale_fill_gradient2(
  low = "#0059FF", high = "#FF1818", mid = "#E8E1DA",
  midpoint = 0, oob = scales::squish, name = "logFC"
)

# --------------------------- base data ----------------------------------
# (1) ALL cells = uncapped; we will use THIS to always color all focal-SN cells
emb_all <- as.data.frame(SingleCellExperiment::reducedDim(milo_obj, "UMAP")) %>%
  setNames(c("UMAP_1","UMAP_2")) %>%
  mutate(cell = colnames(milo_obj)) %>%
  filter(is.finite(UMAP_1), is.finite(UMAP_2))

# (2) background cells = capped version for speed
emb_full <- emb_all
if (nrow(emb_full) > cells_bg_max) {
  set.seed(1)
  emb_full <- emb_full[sample.int(nrow(emb_full), cells_bg_max), ]
}

# Neighborhood UMAP centroids
cent_umap <- cent_umap %>%
  mutate(Nhood = as.integer(Nhood)) %>%
  select(Nhood, UMAP_1, UMAP_2)

# Neighborhood logFC
da_keep <- da %>% transmute(Nhood = as.integer(Nhood), logFC = as.numeric(logFC))

# Membership table (ensure types)
tab_both <- tab_both %>%
  mutate(super_id = as.character(super_id), Nhood = as.integer(Nhood), set = as.character(set)) %>%
  left_join(da_keep, by = "Nhood", suffix = c("", ".da")) %>%
  mutate(logFC = ifelse(is.na(logFC), logFC.da, logFC)) %>% select(-logFC.da)

# Super-centers formatted
super_centers_umap <- super_centers_umap %>% mutate(super_id = as.character(super_id))

# Neighborhood incidence
nh_mat <- milo_obj@nhoods
stopifnot(!is.null(nh_mat), nrow(nh_mat) == ncol(milo_obj))
all_cells <- colnames(milo_obj)

# Helpers to map Nhood IDs -> columns
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

# Precompute each cell’s member neighborhoods (matrix indices)
S <- Matrix::summary(nh_mat)              # (i=row/cell, j=col/nhood) where != 0
cell2nh <- split(S$j, S$i)                # list: cell_row_idx -> integer nhood indices
names(cell2nh) <- as.character(names(cell2nh))

# Map from matrix nhood column index -> Nhood ID as integer (from colnames or position)
nh_colnames <- colnames(nh_mat)
nh_col_to_id <- function(j) as.integer(ifelse(!is.na(suppressWarnings(as.integer(nh_colnames[j]))),
                                              nh_colnames[j], j))

# Quick map: cell -> UMAP coords (uncapped!)
cell_xy <- emb_all %>% select(cell, UMAP_1, UMAP_2)
rownames(cell_xy) <- cell_xy$cell

# --------------------------- core renderer ------------------------------
render_sn_panel <- function(sid, out_dir) {
  stopifnot(sid %in% super_centers_umap$super_id)
  sid_chr <- as.character(sid)
  
  # Neighborhoods in this superneighborhood
  nh_in_sn <- tab_both %>% filter(super_id == sid_chr) %>% pull(Nhood) %>% unique() %>% sort()
  if (!length(nh_in_sn)) return(invisible(NULL))
  
  # Cells belonging to union of those neighborhoods (this is the FULL, uncapped set)
  cells_in_sn <- cells_for_nhoods(nh_in_sn, nh_mat)
  
  # background layers
  emb_bg  <- emb_full
  neigh_bg <- cent_umap %>% filter(is.finite(UMAP_1), is.finite(UMAP_2))
  
  # focal SN neighborhoods (with logFC)
  neigh_sn <- cent_umap %>%
    filter(Nhood %in% nh_in_sn) %>%
    left_join(da_keep, by = "Nhood")
  
  # Super-centroid
  super_xy <- super_centers_umap %>%
    filter(super_id == sid_chr) %>%
    select(super_id, set, super_lfc, UMAP_1, UMAP_2)
  stopifnot(nrow(super_xy) == 1L)
  
  set_i  <- as.character(super_xy$set[1])
  set_col <- pal_set[[set_i]] %||% "#333333"
  
  # Neighborhood → super lines
  seg_neigh2super <- NULL
  if (isTRUE(show_neigh_lines)) {
    seg_neigh2super <- neigh_sn %>%
      transmute(x = UMAP_1, y = UMAP_2,
                xend = super_xy$UMAP_1[1], yend = super_xy$UMAP_2[1])
  }
  
  # Cell → neighborhood (nearest member neighborhood only), with per-neigh + global caps
  seg_cell2neigh <- NULL
  cells_highlight <- NULL
  if (isTRUE(show_cell_lines) || isTRUE(show_cell_points)) {
    nh_xy <- cent_umap %>% filter(Nhood %in% nh_in_sn) %>% select(Nhood, UMAP_1, UMAP_2)
    rownames(nh_xy) <- as.character(nh_xy$Nhood)
    NH <- as.matrix(nh_xy[, c("UMAP_1","UMAP_2"), drop = FALSE])
    
    nearest_member_neigh <- function(cell_idx) {
      j_all <- cell2nh[[as.character(cell_idx)]]
      if (is.null(j_all) || !length(j_all)) return(NA_integer_)
      nh_ids <- nh_col_to_id(j_all)
      nh_ids <- intersect(nh_ids, nh_in_sn)
      if (!length(nh_ids)) return(NA_integer_)
      cbc <- all_cells[cell_idx]
      if (!cbc %in% rownames(cell_xy)) return(NA_integer_)
      cx <- cell_xy[cbc, "UMAP_1"]; cy <- cell_xy[cbc, "UMAP_2"]
      XY <- NH[as.character(nh_ids), , drop = FALSE]
      d  <- sqrt((XY[,1] - cx)^2 + (XY[,2] - cy)^2)
      nh_ids[ which.min(d) ]
    }
    
    cell_idx_all <- match(cells_in_sn, all_cells)
    cell_idx_all <- cell_idx_all[is.finite(cell_idx_all)]
    if (length(cell_idx_all)) {
      assign_df <- tibble::tibble(
        cell_idx = cell_idx_all,
        cell     = all_cells[cell_idx_all]
      ) %>%
        mutate(Nhood = vapply(cell_idx, nearest_member_neigh, integer(1))) %>%
        filter(is.finite(Nhood))
      
      # per-neighborhood cap
      set.seed(1)
      assign_df <- assign_df %>%
        group_by(Nhood) %>%
        mutate(.rand = runif(n())) %>%
        arrange(.rand, .by_group = TRUE) %>%
        slice_head(n = cells_per_neigh_max) %>%
        ungroup() %>%
        select(-.rand)
      
      # global cap
      if (nrow(assign_df) > cells_total_cap) {
        set.seed(1)
        assign_df <- assign_df %>% slice_sample(n = cells_total_cap)
      }
      
      cx <- cell_xy[assign_df$cell, c("UMAP_1","UMAP_2")]
      nx <- nh_xy[as.character(assign_df$Nhood), c("UMAP_1","UMAP_2")]
      seg_cell2neigh <- tibble::tibble(
        x = cx$UMAP_1, y = cx$UMAP_2,
        xend = nx$UMAP_1, yend = nx$UMAP_2
      )
      
      # these are the *sampled* cells we may want to show as bigger points
      cells_highlight <- emb_all %>% semi_join(assign_df, by = "cell")
    }
  }
  
  # ---- NEW: all cells from this SN, uncapped ----
  cells_sn_full <- NULL
  if (length(cells_in_sn)) {
    cells_sn_full <- emb_all %>% filter(cell %in% cells_in_sn)
  }
  
  cap <- glue("{sid_chr}  •  set: {set_i}  •  neighborhoods: {length(nh_in_sn)}  •  cells: {length(cells_in_sn)}")
  
  p <- ggplot() +
    # background cells (capped)
    ggrastr::geom_point_rast(
      data = emb_bg,
      aes(UMAP_1, UMAP_2),
      size = point_sizes$cell_bg, alpha = alphas$cell_bg, color = "grey50", raster.dpi = 600
    ) +
    # background neighborhoods (faint)
    geom_point(
      data = neigh_bg %>% anti_join(neigh_sn, by = "Nhood"),
      aes(UMAP_1, UMAP_2),
      color = "grey60", size = point_sizes$neigh_bg, alpha = alphas$neigh_bg
    ) +
    # cell → neighborhood lines
    { if (!is.null(seg_cell2neigh) && isTRUE(show_cell_lines))
      geom_segment(data = seg_cell2neigh,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   linewidth = 0.2, alpha = alphas$cell_line, color = set_col)
      else NULL } +
    # neighborhood → super lines
    { if (!is.null(seg_neigh2super) && isTRUE(show_neigh_lines))
      geom_segment(data = seg_neigh2super,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   linewidth = 0.5, alpha = alphas$neigh_line, color = set_col)
      else NULL } +
    # focal neighborhoods (filled by logFC)
    geom_point(
      data = neigh_sn,
      aes(UMAP_1, UMAP_2, fill = logFC),
      shape = 21, size = point_sizes$neigh_sn, color = "black",
      alpha = alphas$neigh_sn, stroke = strokes$neigh_sn
    ) +
    bwr_fill +
    # ---- NEW: draw *all* cells of this SN, uncapped, on top of background ----
  { if (!is.null(cells_sn_full) && nrow(cells_sn_full) > 0)
    ggrastr::geom_point_rast(
      data = cells_sn_full,
      aes(UMAP_1, UMAP_2),
      size = point_sizes$cell_sn,
      alpha = alphas$cell_sn,
      color = set_col,
      raster.dpi = 600
    )
    else NULL } +
    # super-centroid (filled by super_lfc, outlined by set color)
    geom_point(
      data = super_xy,
      aes(UMAP_1, UMAP_2, fill = super_lfc),
      shape = 21, size = point_sizes$super, color = set_col,
      alpha = alphas$super, stroke = strokes$super
    ) +
    bwr_fill +
    # label
    ggrepel::geom_text_repel(
      data = super_xy,
      aes(UMAP_1, UMAP_2, label = super_id),
      size = 3.2, seed = 1, min.segment.length = 0.1,
      box.padding = 0.15, max.overlaps = Inf
    ) +
    coord_fixed() +
    theme_void(base_size = 11) +
    theme(legend.position = "top",
          plot.caption = element_text(hjust = 0.5, face = "bold")) +
    labs(caption = cap, fill = "logFC")
  
  pdf_file <- file.path(out_dir, glue("SN_{sid_chr}.pdf"))
  png_file <- file.path(out_dir, glue("SN_{sid_chr}.png"))
  ggsave(pdf_file, p, width = 8.5, height = 8.5)
  ggsave(png_file, p, width = 8.5, height = 8.5, dpi = 300)
  message("Saved: ", pdf_file, " and ", png_file)
  invisible(p)
}

# --------------------------- render all SNs ------------------------------
all_sid <- super_centers_umap$super_id %>% unique()
invisible(lapply(all_sid, render_sn_panel, out_dir = out_dir))

# Combined multi-page PDF
pdf(file.path(out_dir, "ALL_superneighborhoods_2D_layers.pdf"), width = 8.5, height = 8.5)
for (sid in all_sid) print(render_sn_panel(sid, out_dir = out_dir))
dev.off()

message("Done. Outputs in: ", out_dir)
# ======= PLOT: Ridgeline plot of neighborhood logFC ====
suppressPackageStartupMessages({
  library(ggridges)
})

# selection cutoffs
pos_cut <- sig_da %>%
  dplyr::arrange(dplyr::desc(logFC)) %>%
  dplyr::slice(n_top_each) %>%
  dplyr::pull(logFC)

neg_cut <- sig_da %>%
  dplyr::arrange(logFC) %>%
  dplyr::slice(n_top_each) %>%
  dplyr::pull(logFC)

ridge_df <- da2 %>%
  dplyr::filter(is.finite(logFC))

xmax <- max(abs(ridge_df$logFC), na.rm = TRUE)
xlim_use <- c(-xmax, xmax)

col_HC  <- "#0D47A1"  # non-AD-like
col_AD  <- "#C21807"  # AD-like
col_mid <- "grey95"

p_ridge <- ggplot(ridge_df, aes(x = logFC, y = "Neighborhoods")) +
  ggridges::geom_density_ridges_gradient(
    aes(fill = ..x..),
    rel_min_height = 0.001,
    scale = 1.5,
    color = "grey30",
    size  = 0.3
  ) +
  scale_fill_gradient2(
    low      = col_HC,
    mid      = col_mid,
    high     = col_AD,
    midpoint = 0,
    name     = "Neighborhood logFC\n(AD vs non-AD)"
  ) +
  geom_vline(xintercept = neg_cut, linetype = "dashed", linewidth = 0.4, color = col_HC) +
  geom_vline(xintercept = pos_cut, linetype = "dashed", linewidth = 0.4, color = col_AD) +
  coord_cartesian(xlim = xlim_use, expand = TRUE) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    x = "Neighborhood logFC (AD vs non-AD)",
    title = "Non-AD-like → AD-like gradient of neighborhood logFC"
  )

p_ridge
ggsave("Figures/Superneighborhoods/logFC_ridge_top1000_gradient.pdf",
       p_ridge, width = 6, height = 3)

# ======= SAVE BASIC TABLES ==========================================
# SN membership & PCA super-centroids (useful for downstream QC)
readr::write_csv(
  tab_both %>% dplyr::arrange(set, super_id, Nhood),
  "Figures/Superneighborhoods/SN_membership.csv"
)
readr::write_csv(
  super_centers_pca %>% dplyr::arrange(super_id),
  "Figures/Superneighborhoods/SN_super_centroids_PCA.csv"
)

message("Done. Key outputs in Figures/Superneighborhoods/")



# ======= Getting UMAP and PCA coordinates for neighborhoods ==========================================
# get UMAP and PC coordinates for each neighborhood
neigh_coords <- cent_pca %>%
  inner_join(cent_umap, by = "Nhood") %>%           # adds UMAP_1, UMAP_2
  left_join(da %>% transmute(Nhood = as.integer(Nhood), logFC, FDR),
            by = "Nhood")
