## ===== General setup =====
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(purrr)
library(uwot)
library(RColorBrewer)
library(dplyr)

celltypes <- c(
  "Ast.CHI3L1",
  "Ast.DPP10",
  "Ast.GRM3",
  "Exc.L2.3.CBLN2.LINC02306",
  "Exc.L3.4.RORB.CUX2",
  "Exc.L3.5.RORB.PLCH1",
  "Exc.L4.5.RORB.GABRG1",
  "Exc.L4.5.RORB.IL1RAPL2",
  "Exc.L5.6.IT.Car3",
  "Exc.L5.6.NP",
  "Exc.L5.6.RORB.LINC02196",
  "Exc.L5.ET",
  "Exc.L6.CT",
  "Exc.L6.THEMIS.NFIA",
  "Exc.L6b",
  "Exc.NRGN",
  "Exc.RELN.CHD7",
  "CAMs",
  "Mic.P2RY12",
  "Mic.TPT1",
  "T.cells",
  "Inh.ALCAM.TRPM3",
  "Inh.CUX2.MSR1",
  "Inh.ENOX2.SPHKAP",
  "Inh.FBN2.EPB41L4A",
  "Inh.GPC5.RIT2",
  "Inh.L1.2.PAX6.SCGN",
  "Inh.L1.6.LAMP5.CA13",
  "Inh.L1.PAX6.CA4",
  "Inh.L3.5.SST.MAFB",
  "Inh.L5.6.PVALB.STON2",
  "Inh.L5.6.SST.TH",
  "Inh.L6.SST.NPY",
  "Inh.LAMP5.NRG1..Rosehip.",
  "Inh.LAMP5.RELN",
  "Inh.PTPRK.FAM19A1",
  "Inh.PVALB.CA8..Chandelier.",
  "Inh.PVALB.HTR4",
  "Inh.PVALB.SULF1",
  "Inh.RYR3.TSHZ2",
  "Inh.SGCD.PDE3A",
  "Inh.SORCS1.TTN",
  "Inh.VIP.ABI3BP",
  "Inh.VIP.CLSTN2",
  "Inh.VIP.THSD7B",
  "Inh.VIP.TSHZ2",
  "OPCs",
  "Oli",
  "End",
  "Fib.FLRT2",
  "Per",
  "SMC"
)

displayCells <- c( # select cell types for easy of display
  "Inh.L5.6.SST.TH",
  "Inh.ENOX2.SPHKAP",
  "Exc.L3.4.RORB.CUX2",
  "Per",
  "Oli",
  "Ast.DPP10"
)

metaVars <- c(
  "cogn_decline",
  "age_death",
  "msex",
  "amyloid",
  "plaq_n",
  "tangles",
  "gpath"
)

source("./000 utils.r")
meta <- load_meta_data()  # Data frame with projid (donor labels) and AD_status
cell_scores_list <- load_celltype_revised_phenotype("AD_status") # for CPP AD_status score

AD_cellscores <- cell_scores_list[[1]]
AD_cellscoresMat <- as.matrix(AD_cellscores)
pca_result <- prcomp(AD_cellscoresMat, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$projid <- rownames(AD_cellscoresMat)

# merging pca, metadata, and cell scores into one dataframe
merged_df <- merge(meta, pca_df, by.x = "projid", by.y = "projid")
AD_cellscoresWithIds <- AD_cellscores
AD_cellscoresWithIds$projid <- rownames(AD_cellscores)
merged_df_final <- merge(merged_df, AD_cellscoresWithIds, by.x = "projid", by.y = "projid")

# load in inferred AD score (from 005 and Fig 2B)
inferred_ad_phenotype_df<-readRDS('/Users/timschubert/Documents_on_mac/R Projects/MIT/CPP_subtypes/Data/inferred_AD_phenotype.rds')
inferred_ad_phenotype_df$projid <- as.integer(inferred_ad_phenotype_df$projid)
merged_df_final <- merged_df_final %>%
  dplyr::left_join(inferred_ad_phenotype_df, by = "projid")

## ===== Create main heatmap ====

make_donor_heatmap_with_sorting <- function(
    AD_cellscoresMat,
    merged_df_final,
    subgroup        = c("all","non-AD","AD"),
    col_method      = c("kmeans","hclust","none"),
    K               = 4,
    meta_sort_var   = c("tangles","gpath","cogn_decline","AD_status","inferred_AD_phenotype"),
    save_pdf        = TRUE,
    save_mapping    = TRUE,
    out_dir_pdf     = "Results/Clustered_heatmaps/pdf",
    out_dir_csv     = "Results/Clustered_heatmaps/csv",
    out_dir_rds     = "Results/Clustered_heatmaps/rds",
    seed            = 123,
    cluster_colors  = NULL,
    row_mode        = c("cluster","fixed"),
    fixed_row_order = NULL
){
  suppressPackageStartupMessages({
    library(dplyr); library(readr); library(RColorBrewer)
    library(ComplexHeatmap); library(circlize); library(grid)
  })
  subgroup      <- match.arg(subgroup)
  col_method    <- match.arg(col_method)
  meta_sort_var <- match.arg(meta_sort_var)
  row_mode      <- match.arg(row_mode)
  stopifnot(is.matrix(AD_cellscoresMat))
  if (col_method != "hclust") stopifnot(K >= 2)
  
  dir.create(out_dir_pdf, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir_csv, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir_rds, showWarnings = FALSE, recursive = TRUE)
  
  # --- matrix: celltype × donor (display matrix) ---
  mat <- t(AD_cellscoresMat); storage.mode(mat) <- "numeric"
  
  # Subgroup filter
  keep <- switch(
    subgroup,
    "all"    = colnames(mat),
    "AD"     = intersect(colnames(mat), merged_df_final$projid[merged_df_final$AD_status == 1]),
    "non-AD" = intersect(colnames(mat), merged_df_final$projid[merged_df_final$AD_status == 0])
  )
  mat <- mat[, keep, drop = FALSE]
  
  if (nrow(mat) < 2) stop("Too few cell types (rows).")
  if (col_method == "kmeans" && ncol(mat) < K) stop("Not enough donors for requested K in kmeans.")
  if (col_method == "hclust" && ncol(mat) < 2) stop("Need ≥ 2 donors for hclust.")
  
  ann_cols <- c("AD_status","cogn_decline","CR","CDR","amyloid","tangles","gpath")
  ann_base <- merged_df_final[match(colnames(mat), merged_df_final$projid),
                              c("projid", ann_cols), drop = FALSE]
  rownames(ann_base) <- ann_base$projid
  ann_df <- ann_base[, ann_cols, drop = FALSE]
  
  ann_df$AD_status <- factor(ifelse(ann_df$AD_status == 1, "1", "0"), levels = c("0","1"))
  
  # ---------- helpers ----------
  .safe_cor_dist_cols <- function(x) {
    if (ncol(x) <= 1) return(stats::as.dist(matrix(0, nrow = ncol(x), ncol = ncol(x))))
    cmat <- suppressWarnings(stats::cor(x, use = "pairwise.complete.obs", method = "pearson"))
    cmat <- as.matrix(cmat); diag(cmat) <- 1; cmat[!is.finite(cmat)] <- 0
    d <- 1 - cmat; d[!is.finite(d)] <- 1; d <- (d + t(d))/2; d[d < 0] <- 0
    stats::as.dist(d)
  }
  .safe_cor_dist_rows <- function(x) {
    if (nrow(x) <= 1) return(stats::as.dist(matrix(0, nrow = nrow(x), ncol = nrow(x))))
    cmat <- suppressWarnings(stats::cor(t(x), use = "pairwise.complete.obs", method = "pearson"))
    cmat <- as.matrix(cmat); diag(cmat) <- 1; cmat[!is.finite(cmat)] <- 0
    d <- 1 - cmat; d[!is.finite(d)] <- 1; d <- (d + t(d))/2; d[d < 0] <- 0
    stats::as.dist(d)
  }
  .scale_rows <- function(x) { xs <- t(scale(t(x))); xs[!is.finite(xs)] <- 0; xs }
  .apply_fixed_row_order <- function(x, desired) {
    if (is.null(desired)) return(x)
    desired <- unique(desired)
    present <- intersect(desired, rownames(x))
    missing <- setdiff(desired, rownames(x))
    extra   <- setdiff(rownames(x), desired)
    if (length(missing) > 0) message("ℹ Some requested cell types not present: ", paste(missing, collapse = ", "))
    x[c(present, sort(extra)), , drop = FALSE]
  }
  
  # --- cluster assignment (or none) ---
  set.seed(seed)
  kcol <- paste0("K", K)
  cluster_map <- NULL
  original_cluster <- NULL
  
  if (col_method == "kmeans") {
    Xs <- .scale_rows(mat)
    km <- kmeans(t(Xs), centers = K, nstart = 50)
    original_cluster <- km$cluster
    ann_df[[kcol]] <- factor(original_cluster)
  } else if (col_method == "hclust") {
    dcols  <- .safe_cor_dist_cols(mat)
    hc     <- hclust(dcols, method = "ward.D2")
    original_cluster <- cutree(hc, k = K)
    ann_df[[kcol]] <- factor(original_cluster)
  }
  
  # --- sort donors by meta_sort_var ---
  sort_values <- switch(
    meta_sort_var,
    "AD_status"    = ifelse(ann_df$AD_status == "1", 1, 0),
    "tangles"      = suppressWarnings(as.numeric(ann_base$tangles[match(rownames(ann_df), ann_base$projid)])),
    "gpath"        = suppressWarnings(as.numeric(ann_base$gpath[match(rownames(ann_df), ann_base$projid)])),
    "cogn_decline" = suppressWarnings(as.numeric(ann_base$cogn_decline[match(rownames(ann_df), ann_base$projid)]))
  )
  names(sort_values) <- rownames(ann_df)
  
  vis_labels <- NULL
  if (col_method == "none") {
    donors <- rownames(ann_df)
    vals   <- sort_values[donors]
    ord_idx <- order(is.na(vals), vals, donors, na.last = TRUE)
    gaps_col <- NULL
  } else {
    cluster_levels <- levels(ann_df[[kcol]])
    cluster_medians <- sapply(cluster_levels, function(cl){
      vals <- sort_values[ann_df[[kcol]] == cl]
      med  <- suppressWarnings(stats::median(vals, na.rm = TRUE))
      if (!is.finite(med)) Inf else med
    })
    ordered_clusters <- cluster_levels[order(cluster_medians, na.last = TRUE)]
    ord_idx <- unlist(lapply(ordered_clusters, function(cl){
      idx <- which(ann_df[[kcol]] == cl)
      donors <- rownames(ann_df)[idx]
      vals   <- sort_values[donors]
      idx[order(is.na(vals), vals, donors, na.last = TRUE)]
    }), use.names = FALSE)
    
    ann_df[[kcol]] <- factor(as.character(ann_df[[kcol]]), levels = ordered_clusters)
    vis_labels <- paste0("K", seq_along(ordered_clusters))
    levels(ann_df[[kcol]]) <- vis_labels
    
    cl_seq <- as.character(ann_df[[kcol]][ord_idx])
    runlen <- rle(cl_seq)$lengths
    gaps_col <- if (length(runlen) > 1) cumsum(runlen)[-length(runlen)] else NULL
  }
  
  # Reorder matrix & annots by donor order
  mat    <- mat[, ord_idx, drop = FALSE]
  ann_df <- ann_df[ord_idx, , drop = FALSE]
  
  # --- row clustering OR fixed order ---
  if (row_mode == "fixed") {
    mat <- .apply_fixed_row_order(mat, fixed_row_order)
    hr  <- FALSE
  } else {
    drows <- .safe_cor_dist_rows(mat)
    hr    <- hclust(drows, method = "ward.D2")
  }
  
  # ===================== COLORS =====================
  # 1) CPP score colors
  cpp_breaks <- seq(-0.75, 0.75, length.out = 101)
  cpp_cols   <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(101)
  cpp_fun    <- circlize::colorRamp2(cpp_breaks, cpp_cols)
  
  # 2) Annotation colors
  bin_cols <- list(
    AD_status = c("0" = "blue", "1" = "red")
  )
  cont_funs <- list(
    amyloid      = colorRamp2(c(0,   4),   c("white", "#d4a200")),
    tangles      = colorRamp2(c(0,   5),   c("white", "orange")),
    gpath        = colorRamp2(c(0,   4),   c("white", "darkgreen")),
    cogn_decline = colorRamp2(c(-2,  5),   c("white", "purple"))
  )
  
  .build_cluster_colors <- function(cluster_levels, user_cols = NULL) {
    cluster_levels <- as.character(cluster_levels)
    if (!is.null(user_cols)) {
      stopifnot(is.character(user_cols))
      if (!is.null(names(user_cols)) && !all(names(user_cols) %in% c("", NA))) {
        nm <- names(user_cols); nm <- ifelse(grepl("^K\\d+$", nm), nm, paste0("K", nm))
        names(user_cols) <- nm
      } else {
        names(user_cols) <- cluster_levels[seq_len(min(length(user_cols), length(cluster_levels)))]
      }
    } else user_cols <- character(0)
    col_map <- rep(NA_character_, length(cluster_levels)); names(col_map) <- cluster_levels
    overlap <- intersect(names(user_cols), cluster_levels)
    col_map[overlap] <- unname(user_cols[overlap])
    need_n <- sum(is.na(col_map))
    if (need_n > 0) {
      base <- brewer.pal(12,"Set3")
      filler <- if (need_n <= length(base)) base[seq_len(need_n)] else colorRampPalette(base)(need_n)
      col_map[is.na(col_map)] <- filler
    }
    col_map
  }
  if (col_method != "none") {
    bin_cols[[kcol]] <- .build_cluster_colors(levels(ann_df[[kcol]]), user_cols = cluster_colors)
  }
  
  # 3) Top annotation with K* ON TOP, then fixed metadata order
  base_order <- c("AD_status","cogn_decline","CR","CDR","amyloid","tangles","gpath")
  ann_order  <- if (col_method == "none") base_order else c(kcol, base_order)
  
  n_ann <- length(ann_order)
  ha_top <- HeatmapAnnotation(
    df  = ann_df[, ann_order, drop = FALSE],
    col = c(bin_cols, cont_funs),
    simple_anno_size   = unit(1.7, "mm"),
    annotation_height  = rep(unit(1.7, "mm"), n_ann),
    annotation_name_gp = gpar(fontsize = 6),
    gap                = unit(0.4, "mm")
  )
  
  # Column split for visual cluster blocks (adds gaps)
  column_split <- if (col_method == "none") NULL else factor(ann_df[[kcol]], levels = levels(ann_df[[kcol]]))
  
  # --- filenames ---
  cohort_tag <- switch(subgroup, "all" = "ALL", "AD" = "ADONLY", "non-AD" = "HEALTHYONLY")
  method_tag <- switch(col_method,
                       "kmeans" = sprintf("kmeansK%d_sort-%s", K, meta_sort_var),
                       "hclust" = sprintf("hclustK%d_sort-%s", K, meta_sort_var),
                       "none"   = sprintf("nosortKNA_sort-%s",   meta_sort_var))
  pdf_path <- file.path(out_dir_pdf, sprintf("4a.%s.%s.pdf", cohort_tag, method_tag))
  csv_path <- file.path(out_dir_csv, sprintf("4a_donor_clusters.%s.%s.csv", cohort_tag, method_tag))
  rds_path <- file.path(out_dir_rds, sprintf("donor_ids_by_cluster.%s.%s.rds", cohort_tag, method_tag))
  
  # --- draw heatmap ---
  ht <- Heatmap(
    mat,
    name = "CPP",
    col  = cpp_fun,
    cluster_rows = hr,                # FALSE if fixed; hclust if cluster
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp        = gpar(fontsize = 7),
    row_names_max_width = unit(40, "mm"),
    row_title = NULL,
    column_title = sprintf("%s donors • CPP scores • %s (sorted by %s)",
                           switch(subgroup, "all"="All", "AD"="AD-only", "non-AD"="Healthy-only"),
                           gsub("_", " ", method_tag), meta_sort_var),
    top_annotation = ha_top,
    column_split   = column_split
  )
  
  if (isTRUE(save_pdf)) {
    pdf(pdf_path, height = 5, width = 16)
    draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
    dev.off()
    message("✅ Heatmap: ", normalizePath(pdf_path))
  }
  
  # --- save mapping ---
  if (col_method != "none") {
    cluster_map <- data.frame(
      projid               = colnames(mat),
      original_cluster_id  = as.integer(original_cluster[match(colnames(mat), names(original_cluster))]),
      cluster_visual_id    = as.integer(ann_df[[kcol]]),
      cluster_label        = as.character(ann_df[[kcol]]),
      stringsAsFactors     = FALSE
    )
    readr::write_csv(cluster_map, csv_path)
    if (isTRUE(save_mapping)) {
      donors_by <- split(cluster_map$projid, cluster_map$cluster_label)
      saveRDS(donors_by, rds_path)
      message("✅ Mapping CSV/RDS saved.")
    } else {
      rds_path <- NULL
      message("✅ Mapping CSV saved (RDS disabled).")
    }
  } else {
    csv_path <- NULL; rds_path <- NULL
    message("ℹ No discrete column clusters to save (col_method = 'none').")
  }
  
  invisible(list(
    pdf = if (isTRUE(save_pdf)) pdf_path else NULL,
    csv = if (col_method != "none") csv_path else NULL,
    rds = if (col_method != "none" && isTRUE(save_mapping)) rds_path else NULL,
    method = col_method, K = if (col_method == "none") NA_integer_ else K,
    subgroup = subgroup, meta_sort_var = meta_sort_var,
    col_order = colnames(mat),
    annotation = ann_df,
    kcol = if (col_method == "none") NA_character_ else paste0("K", K),
    cluster_map = cluster_map,
    cluster_levels_visual = vis_labels
  ))
}

celltypesorder <- get_fixed_celltype_order_for_plots()

res_all<-make_donor_heatmap_with_sorting(AD_cellscoresMat, merged_df_final,
                                         subgroup = "all", col_method = "kmeans", K = 3,meta_sort_var = "AD_status",
                                         cluster_colors = c("#3B6BA5","#C66E28","#5C946E"),row_mode = "fixed",fixed_row_order = celltypesorder)

# res_ad<-make_donor_heatmap_with_sorting(AD_cellscoresMat, merged_df_final,
#                                 subgroup = "AD", col_method = "kmeans", K = 3, meta_sort_var = "cogn_decline",
#                                 cluster_colors = c("#3B6BA5","#C66E28","#5C946E"),row_mode = "fixed",fixed_row_order = celltypesorder)
# 
# res_nonad<-make_donor_heatmap_with_sorting(AD_cellscoresMat, merged_df_final,
#                                 subgroup = "non-AD", col_method = "kmeans", K = 3, meta_sort_var = "cogn_decline",
#                                 cluster_colors = c("#3B6BA5","#C66E28","#5C946E"),row_mode = "fixed",fixed_row_order = celltypesorder)

# Saving cluster assignments

merge_cluster_assignments <- function(res, meta_df) {
  stopifnot(is.data.frame(res$cluster_map))
  stopifnot("projid" %in% names(res$cluster_map))
  stopifnot("projid" %in% names(meta_df))
  
  # Merge cluster assignments with all metadata columns
  merged <- merge(
    res$cluster_map[, c("projid", "cluster_label")],
    meta_df,
    by = "projid",
    all.x = TRUE
  )
  
  return(merged)
}

df_clusters <- merge_cluster_assignments(res_all, merged_df_final)
head(df_clusters)

# merge into merged_df_final
merged_df_final <- merged_df_final %>%
  mutate(projid = as.character(projid)) %>%
  left_join(
    df_clusters %>% mutate(projid = as.character(projid)) %>% select(projid, cluster_label),
    by = "projid"
  ) %>%
  rename(cluster_labels_for_big_heatmap = cluster_label)


merged_df_final <- merged_df_final %>%
  mutate(
    subgroup = case_when(
      AD_status == 0 & cluster_labels_for_big_heatmap == "K1" ~ "PathNonAD-TxNonAD",
      AD_status == 0 & cluster_labels_for_big_heatmap == "K2" ~ "PathNonAD-TxMid",
      AD_status == 0 & cluster_labels_for_big_heatmap == "K3" ~ "PathNonAD-TxAD",
      AD_status == 1 & cluster_labels_for_big_heatmap == "K1" ~ "PathAD-TxNonAD",
      AD_status == 1 & cluster_labels_for_big_heatmap == "K2" ~ "PathAD-TxMid",
      AD_status == 1 & cluster_labels_for_big_heatmap == "K3" ~ "PathAD-TxAD",
      TRUE ~ NA_character_
    )
  )

table(merged_df_final$subgroup)

saveRDS(merged_df_final,"merged_df_final_11_03.rds")

## ===== Boxplots comparing pathology and clinical outcomes =====
library(ggplot2)
library(dplyr)
library(tidyr)

# Factor ordering
df_clusters$AD_status <- factor(df_clusters$AD_status, levels = c("0", "1"))
df_clusters$cluster_label <- factor(df_clusters$cluster_label, levels = c("K1", "K2", "K3"))

# Create combined group variable (cluster x AD_status)
df_clusters <- df_clusters %>%
  mutate(group = interaction(cluster_label, AD_status, sep = "_"))

# Create a "Total" group for overall distribution per cluster
df_overall <- df_clusters %>%
  mutate(group = as.character(cluster_label),
         AD_status = "overall")

# Combine
df_combined <- bind_rows(df_clusters, df_overall)

# Reshape to long format for plotting multiple variables
df_long <- df_combined %>%
  pivot_longer(cols = c("cogn_decline", "CR", "gpath","pred_AD_prob_glmnet"),
               names_to = "variable",
               values_to = "value")

# Set order for group variable 
group_levels <- c("K1_0","K1_1","K1","K2_0","K2_1","K2","K3_0","K3_1","K3")
df_long$group <- factor(df_long$group, levels = group_levels)

# Set colors
fill_colors <- c("0" = "steelblue", "1" = "firebrick", "overall" = "gray60")

# Plot
ggplot(df_long, aes(x = group, y = value, fill = AD_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
  facet_wrap(~ variable, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = fill_colors) +
  labs(x = "Cluster × AD Status (incl. overall)", y = "Value") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )


# Faceted boxplots for several variables from df_clusters
# one-way ANOVA per (variable, subgroup) with Tukey HSD post-hoc tests

plot_cluster_boxpanels_faceted <- function(
    df_clusters,
    vars                  = c("gpath","cogn_decline","CR"),
    palette               = NULL,
    out_dir               = "Results/Clustered_heatmaps/boxplots",
    out_file_stem         = "boxpanels_ANOVA_Tukey_faceted",
    width                 = 12,
    height                = 6,
    base_size             = 10
){
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2)
    library(ggsignif); library(stringr); library(purrr); library(grid)
  })
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------- helpers ----------
  prettify_cluster <- function(k){
    lab <- c(K1 = "non-AD-like Tx", K2 = "intermediate Tx", K3 = "AD-like Tx")
    kk  <- as.character(k)
    out <- unname(lab[kk])
    out[is.na(out)] <- kk[is.na(out)]
    out
  }
  p_to_stars <- function(p){
    dplyr::case_when(
      is.na(p)  ~ "",
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )
  }
  
  # ---------- input checks ----------
  need_cols <- c("projid","cluster_label","AD_status", vars)
  missing   <- setdiff(need_cols, names(df_clusters))
  if (length(missing) > 0) {
    stop("df_clusters is missing columns: ", paste(missing, collapse = ", "))
  }
  
  # ---------- normalize inputs ----------
  dat <- df_clusters %>%
    transmute(
      projid        = as.character(projid),
      cluster       = case_when(
        cluster_label %in% c("K1","K2","K3") ~ as.character(cluster_label),
        str_detect(cluster_label, "^[Kk]?1$|^K?0?1$") ~ "K1",
        str_detect(cluster_label, "^[Kk]?2$|^K?0?2$") ~ "K2",
        str_detect(cluster_label, "^[Kk]?3$|^K?0?3$") ~ "K3",
        TRUE ~ as.character(cluster_label)
      ),
      AD_status     = factor(as.character(AD_status), levels = c("0","1")),
      !!!rlang::syms(vars)
    ) %>%
    mutate(
      cluster            = factor(cluster, levels = c("K1","K2","K3"), ordered = TRUE),
      cluster_pretty_chr = prettify_cluster(cluster),
      cluster_pretty     = factor(cluster_pretty_chr,
                                  levels = c("non-AD-like Tx","intermediate Tx","AD-like Tx"),
                                  ordered = TRUE)
    ) %>%
    select(-cluster_pretty_chr)
  
  # palette
  if (is.null(palette)) {
    palette <- c(K1="#264653", K2="#f38b20", K3="#ad2524")
  } else {
    stopifnot(all(c("K1","K2","K3") %in% names(palette)))
  }
  
  # ---------- build subgroups: ALL / AD / non-AD ----------
  dat_all    <- dat %>% mutate(subgroup = "ALL")
  dat_ad     <- dat %>% filter(AD_status == "1") %>% mutate(subgroup = "AD")
  dat_nonad  <- dat %>% filter(AD_status == "0") %>% mutate(subgroup = "non-AD")
  
  dat_sub <- bind_rows(dat_all, dat_ad, dat_nonad) %>%
    mutate(
      subgroup = factor(subgroup, levels = c("ALL","AD","non-AD"), ordered = TRUE)
    )
  
  # ---------- long format ----------
  dat_long <- dat_sub %>%
    pivot_longer(all_of(vars), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = vars, ordered = TRUE))
  
  # ---------- Tukey annotations per (variable, subgroup) ----------
  tukey_ann <- dat_long %>%
    filter(!is.na(value)) %>%
    group_by(variable, subgroup) %>%
    group_modify(~{
      df <- .x
      n_per <- df %>% count(cluster)
      if (nrow(n_per %>% filter(n >= 2)) < 2) return(tibble())
      fit <- try(suppressWarnings(aov(value ~ cluster, data = df)), silent = TRUE)
      if (inherits(fit, "try-error")) return(tibble())
      tk <- try(suppressWarnings(TukeyHSD(fit, "cluster")), silent = TRUE)
      if (inherits(tk, "try-error")) return(tibble())
      tk_df <- as.data.frame(tk$cluster)
      tk_df$contrast <- rownames(tk_df)
      out <- tk_df %>%
        transmute(
          contrast,
          group1 = stringr::str_extract(contrast, "^K[123]"),
          group2 = stringr::str_extract(contrast, "K[123]$"),
          p.adj  = `p adj`
        ) %>%
        filter(!is.na(group1), !is.na(group2)) %>%
        mutate(stars = p_to_stars(p.adj)) %>%
        filter(stars != "")
      if (nrow(out) == 0) return(tibble())
      y_max   <- max(df$value, na.rm = TRUE)
      y_min   <- min(df$value, na.rm = TRUE)
      y_range <- ifelse(is.finite(y_max - y_min) && (y_max - y_min) > 0, y_max - y_min, 1)
      step    <- 0.07 * y_range
      out %>%
        arrange(p.adj) %>%
        mutate(
          y.position = y_max + (row_number()) * step,
          xmin = prettify_cluster(group1),
          xmax = prettify_cluster(group2)
        )
    }) %>%
    ungroup()
  
  # ---------- plot ----------
  p <- ggplot(dat_long, aes(x = cluster_pretty, y = value, fill = cluster)) +
    geom_boxplot(width = 0.7, alpha = 0.85, colour = "black", outlier.shape = NA) +
    geom_jitter(aes(color = cluster),
                width = 0.08, height = 0, size = 0.9, alpha = 0.4, shape = 16) +
    scale_fill_manual(values = palette, name = "Cluster (Tx profile)") +
    scale_color_manual(values = palette, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    facet_grid(rows = vars(variable), cols = vars(subgroup), scales = "free_y", switch = "y") +
    theme_classic(base_size = base_size) +
    theme(
      strip.placement  = "outside",
      strip.background = element_rect(fill = "grey95", colour = NA),
      panel.spacing    = unit(0.8, "lines"),
      legend.position  = "right",
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.9)),
      axis.ticks.x     = element_line()
    )
  
  if (nrow(tukey_ann) > 0) {
    p <- p + ggsignif::geom_signif(
      data = tukey_ann,
      aes(xmin = xmin, xmax = xmax, annotations = stars, y_position = y.position),
      manual = TRUE,
      tip_length = 0.01,
      textsize = rel(3.5),
      vjust = -0.1,
      inherit.aes = FALSE
    )
  }
  
  # hide x labels for non-bottom rows
  g <- ggplotGrob(p)
  lay_panels <- g$layout[g$layout$name == "panel", c("t","l","b","r","name")]
  if (nrow(lay_panels) > 0) {
    max_row <- max(lay_panels$t)
    idx_axes_b <- which(grepl("^axis-b-\\d+-\\d+$", g$layout$name))
    for (i in idx_axes_b) {
      nm <- g$layout$name[i]
      rc <- as.integer(stringr::str_match(nm, "^axis-b-(\\d+)-(\\d+)$")[,2])
      row_i <- rc[1]
      if (!is.na(row_i) && row_i < max_row) {
        ax_grob <- g$grobs[[i]]
        try({
          if (!is.null(ax_grob$children[[2]])) ax_grob$children[[2]] <- grid::zeroGrob()
          if (!is.null(ax_grob$grobs) && length(ax_grob$grobs) >= 2) ax_grob$grobs[[2]] <- grid::zeroGrob()
        }, silent = TRUE)
        g$grobs[[i]] <- ax_grob
      }
    }
  }
  
  # save PDF
  outfile <- file.path(out_dir, paste0(out_file_stem, ".pdf"))
  grDevices::pdf(outfile, width = width, height = height, useDingbats = FALSE)
  grid::grid.newpage(); grid::grid.draw(g)
  grDevices::dev.off()
  
  invisible(list(
    plot        = p,
    data        = dat_long,
    annotations = tukey_ann,
    outfile     = outfile
  ))
}


cluster_cols <- c(
  "K1" = "#264653",
  "K2" = "#f38b20",
  "K3" = "#ad2524"
)

bx <- plot_cluster_boxpanels_faceted(
  df_clusters,
  vars          = c("gpath","cogn_decline","CR","pred_AD_prob_glmnet"),
  palette       = cluster_cols,
  out_dir       = "Results/Clustered_heatmaps/boxplots",
  out_file_stem = "boxpanels_ANOVA_Tukey_faceted",
  width         = 5,
  height        = 7,
  base_size     = 10
)



## ===== Differential metabolite analysis on heatmap clusters =====
## ---------- 1) Setup & Loading ---------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel); library(readr)
  library(rlang); library(purrr); library(stringr)
})

## ---------- 2) Input paths and objects ----------

rds_metabolites <- "Data/metabolites_fully_processed.rds"
stopifnot(file.exists(rds_metabolites))
metabolites <- readRDS(rds_metabolites)
stopifnot("projid" %in% names(metabolites))
metabolites$projid <- as.character(metabolites$projid)
stopifnot(!anyDuplicated(metabolites$projid))
met_mat <- as.matrix(metabolites[, setdiff(names(metabolites), "projid"), drop = FALSE])
rownames(met_mat) <- metabolites$projid

stopifnot(!is.null(rownames(met_mat)))
merged_df_final <- readRDS("Data/merged_df_final_11_03.rds")
merged_df_final <- merged_df_final %>% mutate(projid = as.character(projid))

## ---------- 3) Analysis options ----------
OUT_DIR <- "Results/Metabolites"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Covariates
COVARS <- c("pmi","study","fixation_interval","msex")

## ---------- 4) Plot helper ----------
make_volcano <- function(df,
                         estimate_col = "Estimate",
                         p_col        = "p_value",
                         label_col    = "metabolite",
                         title        = NULL,
                         label_fdr    = 0.05,
                         colors       = list(pos = "#D55E00",  
                                             neg = "#0072B2", 
                                             ns  = "grey70"),
                         xlim         = c(-1.5, 1.5),
                         ylim         = c(0, 8),
                         point_size   = 1.8,
                         point_alpha  = 0.85) {
  d <- df
  d$est  <- as.numeric(d[[estimate_col]])
  d$p    <- as.numeric(d[[p_col]])
  d$FDR  <- p.adjust(d$p, method = "fdr")
  d$mlog10FDR <- -log10(pmax(d$FDR, .Machine$double.xmin))
  is_sig <- d$FDR < label_fdr
  
  # map to three classes: pos_sig, neg_sig, ns
  d$sig_dir <- "ns"
  d$sig_dir[is_sig & d$est > 0] <- "pos_sig"
  d$sig_dir[is_sig & d$est < 0] <- "neg_sig"
  
  pal <- c("pos_sig" = colors$pos, "neg_sig" = colors$neg, "ns" = colors$ns)
  
  p <- ggplot(d, aes(x = est, y = mlog10FDR)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    geom_vline(xintercept = 0,            linetype = "solid",  linewidth = 0.5, color = "grey40") +
    geom_point(aes(color = sig_dir), alpha = point_alpha, size = point_size,stroke=0) +
    scale_color_manual(values = pal) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title, x = expression(log[2](FC)), y = expression(-log[10](p[adj]))) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0))
  
  # label only significant points
  lab_df <- d[is_sig & !is.na(d[[label_col]]) & d[[label_col]] != "", , drop = FALSE]
  lab_df$clean_label <- gsub("\\.", " ", lab_df[[label_col]])
  
  if (nrow(lab_df)) {
    suppressPackageStartupMessages(require(ggrepel, quietly = TRUE, warn.conflicts = FALSE))
    p <- p + ggrepel::geom_text_repel(
      data = lab_df[order(lab_df$FDR, -abs(lab_df$est)), , drop = FALSE],
      aes(label = clean_label, color = sig_dir),
      size = 3, max.overlaps = 20, box.padding = 0.35,
      point.padding = 0.25, force = 1, segment.alpha = 0.6, segment.size = 0.3, seed = 123
    )
  }
  
  p
}

## ---------- 5) per-metabolite linear model helper ---------- 
run_binary_comparison_met <- function(met_im, meta, group_df,
                                      ref_level  = NULL,
                                      cohort_ids = NULL,
                                      title_tag  = NULL,
                                      label_fdr  = 0.05,
                                      save_pdf   = TRUE,
                                      pdf_file   = NULL,
                                      save_table = TRUE,
                                      table_file = NULL,
                                      out_dir    = "Results/Group_comparisons",
                                      pdf_width  = 6,
                                      pdf_height = 5,
                                      colors_volcano = list(pos = "#D55E00",
                                                            neg = "#0072B2",
                                                            ns  = "grey70")) {
  message("\n[run_binary_comparison_met] Starting: ", title_tag)
  stopifnot(all(c("projid", "group") %in% names(group_df)))
  stopifnot(all(COVARS %in% names(meta)))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Align IDs
  ids <- Reduce(intersect, list(rownames(met_im), meta$projid, as.character(group_df$projid)))
  if (!is.null(cohort_ids)) ids <- intersect(ids, as.character(cohort_ids))
  if (length(ids) < 3) stop("Too few overlapping samples after filtering.")
  
  X   <- met_im[ids, , drop = FALSE]
  md  <- meta[match(ids, meta$projid), , drop = FALSE]
  grp_raw <- as.character(group_df[match(ids, as.character(group_df$projid)), "group", drop = TRUE])
  
  # Groups and baseline
  levs <- unique(grp_raw)
  if (!ref_level %in% levs) stop(sprintf("ref_level '%s' not present in: %s", ref_level, paste(levs, collapse = ", ")))
  alt_levels <- setdiff(levs, ref_level)
  if (length(alt_levels) != 1) stop(sprintf("Expected exactly 1 alternative level, found: %s", paste(alt_levels, collapse = ", ")))
  alt_level <- alt_levels[1]
  grp <- factor(grp_raw, levels = c(ref_level, alt_level))
  
  message("[run_binary_comparison_met] N total (overlap): ", length(ids))
  message("[run_binary_comparison_met] Ref = ", ref_level, " | Alt = ", alt_level)
  message(sprintf("[run_binary_comparison_met] Group counts: %s=%d | %s=%d",
                  ref_level, sum(grp == ref_level, na.rm = TRUE),
                  alt_level, sum(grp == alt_level, na.rm = TRUE)))
  
  # Covariates (types)
  md$study <- factor(md$study)
  md$msex  <- factor(md$msex)
  md$fixation_interval <- suppressWarnings(as.numeric(md$fixation_interval))
  md$pmi <- suppressWarnings(as.numeric(md$pmi))
  
  # LM per metabolite
  lin_formula_txt <- sprintf("y ~ group + %s", paste(COVARS, collapse = " + "))
  message("[run_binary_comparison_met] Model formula: ", lin_formula_txt)
  message("[run_binary_comparison_met] Running ", ncol(X), " metabolites...")
  
  assoc <- lapply(colnames(X), function(met) {
    df <- data.frame(y = X[, met, drop = TRUE], group = grp, md[, COVARS, drop = FALSE])
    fit <- try(lm(as.formula(lin_formula_txt), data = df, na.action = na.exclude), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(c(metabolite = met, Estimate = NA, StdError = NA, t_value = NA, p_value = NA, n_used = NA))
    }
    co <- summary(fit)$coefficients
    rn <- grep("^group", rownames(co), value = TRUE)[1]
    n_used <- tryCatch(stats::nobs(fit), error = function(...) NA_integer_)
    if (length(rn)) {
      c(metabolite = met,
        Estimate   = unname(co[rn, "Estimate"]),
        StdError   = unname(co[rn, "Std. Error"]),
        t_value    = unname(co[rn, "t value"]),
        p_value    = unname(co[rn, "Pr(>|t|)"]),
        n_used     = n_used)
    } else {
      c(metabolite = met, Estimate = NA, StdError = NA, t_value = NA, p_value = NA, n_used = n_used)
    }
  }) |> do.call(what = "rbind") |> as.data.frame(stringsAsFactors = FALSE)
  
  # Cleanup & FDR correction
  assoc$Estimate <- as.numeric(assoc$Estimate)
  assoc$StdError <- as.numeric(assoc$StdError)
  assoc$t_value  <- as.numeric(assoc$t_value)
  assoc$p_value  <- as.numeric(assoc$p_value)
  assoc$n_used   <- as.integer(assoc$n_used)
  assoc$FDR      <- p.adjust(assoc$p_value, method = "fdr")
  
  n_sig_up <- sum(is.finite(assoc$FDR) & assoc$FDR < label_fdr & assoc$Estimate > 0, na.rm = TRUE)
  message("[run_binary_comparison_met] # significant up (FDR<", label_fdr, "): ", n_sig_up)
  n_sig_down <- sum(is.finite(assoc$FDR) & assoc$FDR < label_fdr & assoc$Estimate < 0, na.rm = TRUE)
  message("[run_binary_comparison_met] # significant down (FDR<", label_fdr, "): ", n_sig_down)
  
  # Plot
  p <- make_volcano(
    assoc,
    estimate_col = "Estimate",
    p_col        = "p_value",
    label_col    = "metabolite",
    title        = sprintf("%s (linear, %s \u2212 %s)", title_tag, alt_level, ref_level),
    label_fdr    = label_fdr,
    colors       = colors_volcano
  )
  
  saved_pdf <- NULL
  if (isTRUE(save_pdf)) {
    if (is.null(pdf_file)) {
      safe <- gsub("[^A-Za-z0-9._-]+", "_", tolower(title_tag))
      pdf_file <- file.path(out_dir, sprintf("%s_volcano.pdf", safe))
    }
    ggsave(pdf_file, plot = p, device = "pdf", width = pdf_width, height = pdf_height, units = "in")
    saved_pdf <- normalizePath(pdf_file)
    message("[run_binary_comparison_met] Saved PDF: ", saved_pdf)
  }
  
  saved_csv <- NULL
  if (isTRUE(save_table)) {
    if (is.null(table_file)) {
      safe <- gsub("[^A-Za-z0-9._-]+", "_", tolower(title_tag))
      table_file <- file.path(out_dir, sprintf("%s_linear_results.csv", safe))
    }
    readr::write_csv(assoc, table_file)
    saved_csv <- normalizePath(table_file)
    message("[run_binary_comparison_met] Saved table: ", saved_csv)
  }
  
  attr(assoc, "contrast_info") <- list(
    ref   = ref_level,
    alt   = alt_level,
    covars = COVARS,
    effect = sprintf("Estimate = mean(log-intensity) in %s minus %s (adjusted)", alt_level, ref_level),
    pdf   = saved_pdf, csv = saved_csv, n = length(ids)
  )
  
  list(assoc_table = assoc, plot = p, pdf = saved_pdf, csv = saved_csv)
}

## ---------- 6) Build groups ---------- 
# a) PathNonAD-TxNon vs. PathAD-TxAD (baseline PathNonAD-TxNon)
group_df_subgroups <- merged_df_final %>%
  transmute(projid = as.character(projid),
            group = as.character(subgroup)) %>%
  filter(!is.na(group) & group %in% c("PathNonAD-TxNonAD", "PathAD-TxAD"))

# b) AD 1 vs 0 (baseline 0)
group_df_AD <- merged_df_final %>%
  transmute(projid = as.character(projid),
            group  = as.character(AD_status)) %>%
  filter(!is.na(group) & group %in% c("0","1"))

# c) K3 vs K1 (baseline K1)
group_df_clusters <- merged_df_final %>%
  transmute(projid = as.character(projid),
            group  = as.character(cluster_labels_for_big_heatmap)) %>%
  filter(!is.na(group) & group %in% c("K1","K3"))

# d) discordant [PathNonAD-TxAD vs. PathAd-TxNon](baseline PathNonAD-TxAD)
group_df_discordant <- merged_df_final %>%
  transmute(projid = as.character(projid),
            group = as.character(subgroup)) %>%
  filter(!is.na(group) & group %in% c("PathAD-TxNonAD", "PathNonAD-TxAD"))

## 6b get number of donors in each group who also have metabolites ##
ids_overlap <- intersect(rownames(met_mat), group_df_subgroups$projid)
group_overlap <- group_df_subgroups$group[match(ids_overlap, group_df_subgroups$projid)]
table(group_overlap)

ids_overlap <- intersect(rownames(met_mat), group_df_AD$projid)
group_overlap <- group_df_AD$group[match(ids_overlap, group_df_AD$projid)]
table(group_overlap)

ids_overlap <- intersect(rownames(met_mat), group_df_clusters$projid)
group_overlap <- group_df_clusters$group[match(ids_overlap, group_df_clusters$projid)]
table(group_overlap)

ids_overlap <- intersect(rownames(met_mat), group_df_discordant$projid)
group_overlap <- group_df_discordant$group[match(ids_overlap, group_df_discordant$projid)]
table(group_overlap)

## ---------- 7) Run contrasts ---------- 

meta_for_covars <- merged_df_final %>%
  select(projid, all_of(COVARS)) %>%
  mutate(projid = as.character(projid))

# Ensure covariates exist
missing_cov <- setdiff(COVARS, names(meta_for_covars))
if (length(missing_cov)) stop("Missing covariates in merged_df_final: ", paste(missing_cov, collapse = ", "))

# Colors
cols_k <- list(pos = "#ad2524", neg = "#264653", ns = "grey80")
cols_g <- list(pos = "#ad2524", neg = "#264653", ns = "grey80")  # PathAD-TxAD, PathNonAD-TxNonAD
cols_a <- list(pos = "#F2696B", neg = "#869AD5", ns = "grey80")  # AD, NonAD

# a) PathNonAD-TxNon vs. PathAD-TxAD (baseline PathNonAD-TxNon)
res_subgroups <- run_binary_comparison_met(
  met_im     = met_mat,
  meta       = meta_for_covars,
  group_df   = group_df_subgroups,
  ref_level  = "PathNonAD-TxNonAD",
  title_tag  = "Subgroups PathNonAD-TxNonAD vs. PathAD-TxAD",
  label_fdr  = 0.05,
  save_pdf   = TRUE,
  pdf_file   = file.path(OUT_DIR, "subgroups_PathNonADTxNon_vs_PathADTxAD_volcano.pdf"),
  save_table = TRUE,
  table_file = file.path(OUT_DIR, "subgroups_PathNonADTxNon_vs_PathADTxAD_results.csv"),
  out_dir    = OUT_DIR,
  pdf_width  = 6, pdf_height = 5,
  colors_volcano = cols_g
)

# b) AD 1 vs 0 (baseline 0)
res_ad_vs_ctrl <- run_binary_comparison_met(
  met_im     = met_mat,
  meta       = meta_for_covars,
  group_df   = group_df_AD,
  ref_level  = "0",
  title_tag  = "AD 1 vs 0",
  label_fdr  = 0.05,
  save_pdf   = TRUE,
  pdf_file   = file.path(OUT_DIR, "AD_1_vs_0_volcano.pdf"),
  save_table = TRUE,
  table_file = file.path(OUT_DIR, "AD_1_vs_0_linear_results.csv"),
  out_dir    = OUT_DIR,
  pdf_width  = 6, pdf_height = 5,
  colors_volcano = cols_a
)

# c) K3 vs K1 (baseline K1)
res_k3_vs_k1 <- run_binary_comparison_met(
  met_im     = met_mat,
  meta       = meta_for_covars,
  group_df   = group_df_clusters,
  ref_level  = "K1",
  title_tag  = "Clusters K3 vs K1",
  label_fdr  = 0.05,
  save_pdf   = TRUE,
  pdf_file   = file.path(OUT_DIR, "clusters_k3_vs_k1_volcano.pdf"),
  save_table = TRUE,
  table_file = file.path(OUT_DIR, "clusters_k3_vs_k1_linear_results.csv"),
  out_dir    = OUT_DIR,
  pdf_width  = 6, pdf_height = 5,
  colors_volcano = cols_k
)

# d) K3 vs K1 among AD cases only (AD_status == 1)
ad1_ids <- merged_df_final %>%
  filter(!is.na(AD_status) & AD_status == 1) %>%
  transmute(projid = as.character(projid)) %>%
  distinct() %>%
  pull(projid)

res_k3_vs_k1_AD1 <- run_binary_comparison_met(
  met_im       = met_mat,
  meta         = meta_for_covars,           
  group_df     = group_df_clusters,         # filtering via cohort_ids below
  ref_level    = "K1",
  cohort_ids   = ad1_ids,                   # restrict to AD==1 donors
  title_tag    = "Clusters K3 vs K1 (AD=1 only)",
  label_fdr    = 0.05,
  save_pdf     = TRUE,
  pdf_file     = file.path(OUT_DIR, "clusters_k3_vs_k1_AD1_volcano.pdf"),
  save_table   = TRUE,
  table_file   = file.path(OUT_DIR, "clusters_k3_vs_k1_AD1_linear_results.csv"),
  out_dir      = OUT_DIR,
  pdf_width    = 6,
  pdf_height   = 5,
  colors_volcano = cols_k                  
)

# e) K3 vs K1 among non-AD cases only 
nonad_ids <- merged_df_final %>%
  dplyr::filter(!is.na(AD_status) & AD_status == 0) %>%
  dplyr::transmute(projid = as.character(projid)) %>%
  dplyr::distinct() %>%
  dplyr::pull(projid)

res_k3_vs_k1_nonAD <- run_binary_comparison_met(
  met_im       = met_mat,
  meta         = meta_for_covars,
  group_df     = group_df_clusters, 
  ref_level    = "K1",
  cohort_ids   = nonad_ids,           # <<< restrict to healthy / AD==0
  title_tag    = "Clusters K3 vs K1 (AD=0 only)",
  label_fdr    = 0.05,
  save_pdf     = TRUE,
  pdf_file     = file.path(OUT_DIR, "clusters_k3_vs_k1_AD0_volcano.pdf"),
  save_table   = TRUE,
  table_file   = file.path(OUT_DIR, "clusters_k3_vs_k1_AD0_linear_results.csv"),
  out_dir      = OUT_DIR,
  pdf_width    = 6,
  pdf_height   = 5,
  colors_volcano = cols_k       
)

# f) discordant [PathNonAD-TxAD vs. PathAd-TxNon](baseline PathNonAD-TxAD)
res_discordant <- run_binary_comparison_met(
  met_im     = met_mat,
  meta       = meta_for_covars,
  group_df   = group_df_discordant,
  ref_level  = "PathNonAD-TxAD",
  title_tag  = "Discordant Subgroups PathNonAD-TxAD vs. PathAd-TxNon",
  label_fdr  = 0.05,
  save_pdf   = TRUE,
  pdf_file   = file.path(OUT_DIR, "discordant_subgroups_PathNonADTxAD_vs_PathAdTxNon_volcano.pdf"),
  save_table = TRUE,
  table_file = file.path(OUT_DIR, "discordant_subgroups_PathNonADTxAD_vs_PathAdTxNon_results.csv"),
  out_dir    = OUT_DIR,
  pdf_width  = 6, pdf_height = 5,
  colors_volcano = cols_g
)

## ---------- 9) Console messages ---------- 
if (!is.null(res_k3_vs_k1$assoc_table)) {
  message("\nTop K3 vs K1 (by FDR):")
  print(res_k3_vs_k1$assoc_table %>% arrange(FDR) %>% slice_head(n = 10))
}
if (!is.null(res_ad_vs_ctrl$assoc_table)) {
  message("\nTop AD 1 vs 0 (by FDR):")
  print(res_ad_vs_ctrl$assoc_table %>% arrange(FDR) %>% slice_head(n = 10))
}
## ---------- 10) Effect size comparison scatter plot ------
plot_two_comparisons_scatter <- function(
    file_x, file_y,
    name_x, name_y,                 
    out_dir        = OUT_DIR,
    fdr_cut        = 0.05,
    alpha_sig      = 0.8,
    alpha_nonsig   = 0.2,
    size_pts       = 1.7,
    seed_labels    = 123,
    pdf_name       = NULL,
    label_mode     = c("smart","sig","all","none"),
    max_labels     = 25,
    force_labels   = NULL,
    lim            = 1.5,            
    clean_labels   = TRUE,
    label_clean_fun = function(x) gsub("\\.", " ", x),
    colors = list(
      both = "black",
      up_x_only = "#D55E00",
      down_x_only = "#E69F00",
      up_y_only = "#0072B2",
      down_y_only = "#56B4E9",
      ns = "grey50"
    )
){
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(ggrepel); library(readr) })
  label_mode <- match.arg(label_mode)
  stopifnot(file.exists(file_x), file.exists(file_y))
  
  message("\n[scatter] Reading:")
  message("  X: ", normalizePath(file_x))
  message("  Y: ", normalizePath(file_y))
  
  tx <- readr::read_csv(file_x, show_col_types = FALSE)
  ty <- readr::read_csv(file_y, show_col_types = FALSE)
  need <- c("metabolite","Estimate","FDR")
  stopifnot(all(need %in% names(tx)), all(need %in% names(ty)))
  
  df <- tx %>%
    select(metabolite, est_x = Estimate, fdr_x = FDR) %>%
    inner_join(ty %>% select(metabolite, est_y = Estimate, fdr_y = FDR), by = "metabolite") %>%
    mutate(sig_x = fdr_x < fdr_cut, sig_y = fdr_y < fdr_cut)
  
  message("[scatter] Overlap metabolites: ", nrow(df))
  message(sprintf("[scatter] #sig (FDR<%.2g): %s=%d, %s=%d",
                  fdr_cut, name_x, sum(df$sig_x, na.rm = TRUE), name_y, sum(df$sig_y, na.rm = TRUE)))
  
  # ---- Significance & direction categories ----
  df <- df %>%
    mutate(
      sig_cat = case_when(
        sig_x & sig_y ~ "both",
        sig_x & !sig_y & est_x > 0 ~ "up_x_only",
        sig_x & !sig_y & est_x < 0 ~ "down_x_only",
        sig_y & !sig_x & est_y > 0 ~ "up_y_only",
        sig_y & !sig_x & est_y < 0 ~ "down_y_only",
        TRUE ~ "ns"
      ),
      alpha = ifelse(sig_cat == "ns", alpha_nonsig, alpha_sig),
      abs_sum = abs(est_x) + abs(est_y),
      quadrant = paste0(ifelse(est_x >= 0, "R", "L"), ifelse(est_y >= 0, "U", "D"))
    )
  
  lvl <- c("ns","up_y_only","down_y_only","up_x_only","down_x_only","both")
  stopifnot(all(lvl %in% names(colors)))
  df$sig_cat <- factor(df$sig_cat, levels = lvl)
  df <- df[order(df$sig_cat), , drop = FALSE]  # background first, foreground last
  
  # ---- Console summary ----
  message("[scatter] Category counts:")
  print(table(df$sig_cat, useNA = "ifany"))
  
  # ---- Label selection ----
  label_idx <- rep(FALSE, nrow(df)); names(label_idx) <- df$metabolite
  if (!is.null(force_labels) && length(force_labels)) {
    label_idx <- label_idx | (df$metabolite %in% force_labels)
  }
  if (label_mode == "all") {
    label_idx[] <- TRUE
  } else if (label_mode == "none") {
    label_idx[] <- FALSE
  } else if (label_mode == "sig") {
    cand <- which(df$sig_cat != "ns")
    # order by our category order, then by effect size
    cand <- cand[order(match(as.character(df$sig_cat[cand]), lvl), -df$abs_sum[cand])]
    label_idx[head(cand, max_labels)] <- TRUE
  } else { 
    idx_sig <- which(df$sig_cat != "ns")
    ord_sig <- idx_sig[order(-df$abs_sum[idx_sig])]
    chosen <- ord_sig
    if (length(chosen) > max_labels) {
      keep <- integer(0)
      for (q in c("RU","LU","RD","LD")) {
        qidx <- chosen[df$quadrant[chosen] == q]
        if (length(qidx)) keep <- c(keep, head(qidx, 2))
      }
      remain <- setdiff(chosen, keep)
      keep <- unique(c(keep, head(remain[order(-df$abs_sum[remain])], max(0, max_labels - length(keep)))))
      chosen <- keep
    }
    label_idx[chosen] <- TRUE
  }
    idx_both <- which(df$sig_cat == "both")
  if (length(idx_both) >= 1) {
    top_both <- idx_both[ which.max(df$abs_sum[idx_both]) ]
    label_idx[top_both] <- TRUE
  }
  
  df$label_plot <- df$metabolite
  if (isTRUE(clean_labels)) {
    df$label_plot <- label_clean_fun(df$label_plot)
  }
  
  # ---- Correlation ----
  ok <- is.finite(df$est_x) & is.finite(df$est_y)
  cor_xy <- if (sum(ok) >= 3) suppressWarnings(stats::cor(df$est_x[ok], df$est_y[ok], method = "pearson")) else NA
  message(sprintf("[scatter] Pearson r = %.3f (n=%d)", cor_xy, sum(ok)))
  
  # ---- Plot ----
  title_txt <- sprintf("%s vs %s", name_x, name_y)
  
  p <- ggplot(df, aes(x = est_x, y = est_y, color = sig_cat, alpha = alpha)) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.4, color = "grey60") +
    geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.4, color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
    geom_point(size = size_pts, stroke = 0) +
    scale_color_manual(
      breaks = lvl,
      values = unlist(colors)[lvl],
      drop = FALSE
    ) +
    scale_alpha_identity() +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = TRUE) +
    ggrepel::geom_text_repel(
      data = df[label_idx, , drop = FALSE],
      aes(label = label_plot),
      size = 3, max.overlaps = 1000, box.padding = 0.4, point.padding = 0.3,
      force = 2.0, min.segment.length = 0, segment.alpha = 0.6, segment.size = 0.3, seed = seed_labels
    ) +
    labs(title = title_txt,
         x = bquote(log[2]*" FC ("*.(name_x)*")"),
         y = bquote(log[2]*" FC ("*.(name_y)*")")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.box = "horizontal", plot.title = element_text(face = "bold", hjust = 0.5))
  
  if (is.null(pdf_name)) {
    base <- paste0(gsub("[^A-Za-z0-9._-]+","_",name_x), "_vs_",
                   gsub("[^A-Za-z0-9._-]+","_",name_y), "_scatter.pdf")
    pdf_name <- file.path(out_dir, base)
  }
  ggsave(pdf_name, plot = p, device = "pdf", width = 5, height = 5, units = "in", dpi = 300)
  message("✅ Saved scatter: ", normalizePath(pdf_name))
  invisible(list(plot = p, data = df, pdf = pdf_name, cor = cor_xy))
}

file_k  <- file.path(OUT_DIR, "clusters_k3_vs_k1_linear_results.csv")
file_c <- file.path(OUT_DIR, "subgroups_PathNonADTxNon_vs_PathADTxAD_results.csv")
file_ad <- file.path(OUT_DIR, "AD_1_vs_0_linear_results.csv")



# a) AD vs Control / K3 vs K1

my_cols <- list(
  both = "black",
  up_x_only = "#F2696B",    
  down_x_only = "#869AD5",  
  up_y_only = "#AD2524",  
  down_y_only = "#264653", 
  ns = "grey50"
)

scatter_comp <- plot_two_comparisons_scatter(
  file_x     = file_ad,
  file_y     = file_k,
  name_x     = "AD (1) \u2212 Ctrl (0)",
  name_y     = "K3 \u2212 K1",
  out_dir    = OUT_DIR,
  fdr_cut    = 0.05,
  label_mode = "smart",
  max_labels = 25,
  lim        = 1.5,
  colors     = my_cols,
  clean_labels = TRUE,
  alpha_nonsig = 0.17,
)

# b) AD vs. Control / PathAD;TxAD vs. PathNonAD;TxNon
subgroup_scatter_comp <- plot_two_comparisons_scatter(
  file_x     = file_ad,
  file_y     = file_c,
  name_x     = "AD (1) \u2212 Ctrl (0)",
  name_y     = "PathAD;TxAD \u2212 PathNonAD;TxNon",
  out_dir    = OUT_DIR,
  fdr_cut    = 0.05,
  label_mode = "smart",
  max_labels = 25,
  lim        = 1.5,
  colors     = my_cols,
  clean_labels = TRUE,
  alpha_nonsig = 0.17,
)

# c) K3 vs K1 among AD cases / K3 vs K1 in all

file_k_ad_pos<-file.path(OUT_DIR,"clusters_k3_vs_k1_AD1_linear_results.csv")

my_cols <- list(
  both = "black",
  up_x_only = "#7B497D",    
  down_x_only = "#C8857F",  
  up_y_only = "#AD2524",  
  down_y_only = "#264653", 
  ns = "grey50"
)

scatter_comp <- plot_two_comparisons_scatter(
  file_x     = file_k_ad_pos,
  file_y     = file_k,
  name_x     = "K3 vs K1 in AD",
  name_y     = "K3 vs K1",
  out_dir    = OUT_DIR,
  fdr_cut    = 0.05,
  label_mode = "smart",
  max_labels = 100,
  lim        = 1.5,
  colors     = my_cols,
  clean_labels = TRUE,
  alpha_nonsig = 0.17,
)

#d) subgroups vs. discordant
file_discordant_pos <- file.path(OUT_DIR, "discordant_subgroups_PathNonADTxAD_vs_PathAdTxNon_results.csv")

discordant_scatter_comp <- plot_two_comparisons_scatter(
  file_x     = file_c,
  file_y     = file_discordant_pos,
  name_x     = "Subgroups",
  name_y     = "Discordant",
  out_dir    = OUT_DIR,
  fdr_cut    = 0.05,
  label_mode = "smart",
  max_labels = 10,
  lim        = 1.5,
  colors     = my_cols,
  clean_labels = TRUE,
  alpha_nonsig = 0.17,
)




















## ===== Cowplot combining scatter + volcanoes (SUBGROUP) ==============
suppressPackageStartupMessages({ library(cowplot); library(ggplot2); library(dplyr); library(ggrepel) })

# Shared constants
LIM_FC  <- 1.5
BRKS_FC <- seq(-1.5, 1.5, by = 0.5)
ALPHA_NS  <- 0.17
ALPHA_SIG <- 0.85

ZERO_LINE_COLOR <- "grey60"; ZERO_LINE_SIZE <- 0.45
DIAG_LINE_COLOR <- "grey55"; DIAG_LINE_SIZE <- 0.55
FDR_LINE_COLOR  <- "grey50"; FDR_LINE_SIZE  <- 0.45
m_tight <- margin(4, 6, 4, 6)

# AD + cluster data 
ad_df <- res_ad_vs_ctrl$assoc_table %>%
  transmute(
    metabolite,
    est = as.numeric(Estimate),
    p   = as.numeric(p_value),
    FDR = p.adjust(p, method = "fdr"),
    mlog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    sig_dir = case_when(FDR < 0.05 & est > 0 ~ "pos_sig",
                        FDR < 0.05 & est < 0 ~ "neg_sig",
                        TRUE ~ "ns"),
    alpha_val  = ifelse(sig_dir == "ns", ALPHA_NS, ALPHA_SIG),
    clean_label = gsub("\\.", " ", metabolite)
  )

subgroup_df <-res_subgroups$assoc_table %>%
  transmute(
    metabolite,
    est = as.numeric(Estimate),
    p   = as.numeric(p_value),
    FDR = p.adjust(p, method = "fdr"),
    mlog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    sig_dir = case_when(FDR < 0.05 & est > 0 ~ "pos_sig",
                        FDR < 0.05 & est < 0 ~ "neg_sig",
                        TRUE ~ "ns"),
    alpha_val  = ifelse(sig_dir == "ns", ALPHA_NS, ALPHA_SIG),
    clean_label = gsub("\\.", " ", metabolite)
  )

# Shared significance limits
SIG_MAX  <- suppressWarnings(max(c(ad_df$mlog10FDR, subgroup_df$mlog10FDR), na.rm = TRUE))
if (!is.finite(SIG_MAX) || SIG_MAX <= 0) SIG_MAX <- 5
SIG_BRKS <- seq(0, ceiling(SIG_MAX), by = 1)
SIG_LIMS <- c(0, ceiling(SIG_MAX))

pal_ad <- c(pos_sig="#b30059", neg_sig="#0066cc", ns="grey60")
pal_subgroup  <- c(pos_sig="#ad2524", neg_sig="#264653", ns="grey60")
pal_sc <- c(both="black", up_x_only="#F2696B", down_x_only="#869AD5",
            up_y_only="#AD2524", down_y_only="#264653", ns="grey60")

# Helper: top3 up/down 
pick_top3_labels <- function(df) {
  up   <- df %>% filter(FDR < 0.05, est > 0) %>% arrange(desc(est)) %>% slice_head(n = 3)
  down <- df %>% filter(FDR < 0.05, est < 0) %>% arrange(est)        %>% slice_head(n = 3)
  bind_rows(up, down)
}
lab_ad <- pick_top3_labels(ad_df)
lab_subgroup  <- pick_top3_labels(subgroup_df)

# AD volcano 
ad_vol <- ggplot(ad_df, aes(x = est, y = mlog10FDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = FDR_LINE_SIZE, color = FDR_LINE_COLOR) +
  geom_vline(xintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_point(aes(color = sig_dir, alpha = alpha_val), size = 1.8, stroke = 0) +
  scale_color_manual(values = pal_ad, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  scale_y_continuous(limits = SIG_LIMS, breaks = SIG_BRKS, expand = expansion(mult = 0)) +
  ggrepel::geom_text_repel(
    data = lab_ad, aes(label = clean_label, color = sig_dir),
    size = 3, box.padding = 0.35, point.padding = 0.25, max.overlaps = 1000,
    force = 1, segment.alpha = 0.6, segment.size = 0.3, seed = 123, show.legend = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    plot.title      = element_blank(),
    plot.margin     = m_tight
  )

# subgroup volcano (rotated left) 
subgroup_vol_rot <- ggplot(subgroup_df, aes(x = est, y = mlog10FDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = FDR_LINE_SIZE, color = FDR_LINE_COLOR) +
  geom_vline(xintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_point(aes(color = sig_dir, alpha = alpha_val), size = 1.8, stroke = 0) +
  scale_color_manual(values = pal_subgroup, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  coord_flip() +
  scale_y_reverse(limits = rev(SIG_LIMS), breaks = rev(SIG_BRKS), expand = expansion(mult = 0)) +
  ggrepel::geom_text_repel(
    data = lab_subgroup, aes(label = clean_label, color = sig_dir),
    size = 3, box.padding = 0.35, point.padding = 0.25, max.overlaps = 1000,
    force = 1, segment.alpha = 0.6, segment.size = 0.3, seed = 123, show.legend = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    plot.title      = element_blank(),
    plot.margin     = m_tight
  )

# Scatter
sc_df <- subgroup_scatter_comp$data %>%
  mutate(alpha_val = ifelse(sig_cat == "ns", ALPHA_NS, ALPHA_SIG))

# labels
label_idx <- rep(FALSE, nrow(sc_df))
idx_sig   <- which(sc_df$sig_cat != "ns")
ord_sig   <- idx_sig[order(-sc_df$abs_sum[idx_sig])]
keep <- integer(0)
for (q in c("RU","LU","RD","LD")) {
  qidx <- ord_sig[sc_df$quadrant[ord_sig] == q]
  if (length(qidx)) keep <- c(keep, head(qidx, 2))
}
remain <- setdiff(ord_sig, keep)
keep <- unique(c(keep, head(remain[order(-sc_df$abs_sum[remain])], max(0, 25 - length(keep)))))
idx_both <- which(sc_df$sig_cat == "both")
if (length(idx_both) >= 1) keep <- unique(c(keep, idx_both[which.max(sc_df$abs_sum[idx_both])]))
label_idx[keep] <- TRUE

scatter_p <- ggplot(sc_df, aes(x = est_x, y = est_y, color = sig_cat, alpha = alpha_val)) +
  geom_hline(yintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_vline(xintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = DIAG_LINE_SIZE, color = DIAG_LINE_COLOR) +
  geom_point(size = 1.7, stroke = 0) +
  scale_color_manual(values = pal_sc, breaks = names(pal_sc), drop = FALSE, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  coord_equal() +
  ggrepel::geom_text_repel(
    data = sc_df[label_idx, , drop = FALSE],
    aes(label = label_plot),
    size = 3, box.padding = 0.4, point.padding = 0.3,
    force = 2.0, min.segment.length = 0, segment.alpha = 0.6, segment.size = 0.3, seed = 123
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    plot.title      = element_blank(),
    plot.margin     = m_tight
  )

# Layout
g_combo <- ggdraw() +
  draw_plot(subgroup_vol_rot, x = 0/3, y = 0/3, width = 1/3, height = 2/3) +
  draw_plot(ad_vol,    x = 1/3, y = 2/3, width = 2/3, height = 1/3) +
  draw_plot(scatter_p, x = 1/3, y = 0/3, width = 2/3, height = 2/3)

# Save
out_pdf <- file.path(OUT_DIR, "AD_subgroup_scatter_cowplot_grid_ticks_top3.pdf")
ggsave(out_pdf, plot = g_combo, device = "pdf", width = 7.5, height = 7.5, units = "in")
message("✅ Saved:", normalizePath(out_pdf))

























## ===== Cowplot combining scatter + volcanoes (CLUSTER) ==============
suppressPackageStartupMessages({ library(cowplot); library(ggplot2); library(dplyr); library(ggrepel) })

# Shared constants
LIM_FC  <- 1.5
BRKS_FC <- seq(-1.5, 1.5, by = 0.5)
ALPHA_NS  <- 0.17
ALPHA_SIG <- 0.85

ZERO_LINE_COLOR <- "grey60"; ZERO_LINE_SIZE <- 0.45
DIAG_LINE_COLOR <- "grey55"; DIAG_LINE_SIZE <- 0.55
FDR_LINE_COLOR  <- "grey50"; FDR_LINE_SIZE  <- 0.45
m_tight <- margin(4, 6, 4, 6)

# AD + cluster data 
ad_df <- res_ad_vs_ctrl$assoc_table %>%
  transmute(
    metabolite,
    est = as.numeric(Estimate),
    p   = as.numeric(p_value),
    FDR = p.adjust(p, method = "fdr"),
    mlog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    sig_dir = case_when(FDR < 0.05 & est > 0 ~ "pos_sig",
                        FDR < 0.05 & est < 0 ~ "neg_sig",
                        TRUE ~ "ns"),
    alpha_val  = ifelse(sig_dir == "ns", ALPHA_NS, ALPHA_SIG),
    clean_label = gsub("\\.", " ", metabolite)
  )

k_df <- res_k3_vs_k1$assoc_table %>%
  transmute(
    metabolite,
    est = as.numeric(Estimate),
    p   = as.numeric(p_value),
    FDR = p.adjust(p, method = "fdr"),
    mlog10FDR = -log10(pmax(FDR, .Machine$double.xmin)),
    sig_dir = case_when(FDR < 0.05 & est > 0 ~ "pos_sig",
                        FDR < 0.05 & est < 0 ~ "neg_sig",
                        TRUE ~ "ns"),
    alpha_val  = ifelse(sig_dir == "ns", ALPHA_NS, ALPHA_SIG),
    clean_label = gsub("\\.", " ", metabolite)
  )

# Shared significance limits
SIG_MAX  <- suppressWarnings(max(c(ad_df$mlog10FDR, k_df$mlog10FDR), na.rm = TRUE))
if (!is.finite(SIG_MAX) || SIG_MAX <= 0) SIG_MAX <- 5
SIG_BRKS <- seq(0, ceiling(SIG_MAX), by = 1)
SIG_LIMS <- c(0, ceiling(SIG_MAX))

pal_ad <- c(pos_sig="#b30059", neg_sig="#0066cc", ns="grey60")
pal_k  <- c(pos_sig="#ad2524", neg_sig="#264653", ns="grey60")
pal_sc <- c(both="black", up_x_only="#F2696B", down_x_only="#869AD5",
            up_y_only="#AD2524", down_y_only="#264653", ns="grey60")

# Helper: top3 up/down 
pick_top3_labels <- function(df) {
  up   <- df %>% filter(FDR < 0.05, est > 0) %>% arrange(desc(est)) %>% slice_head(n = 3)
  down <- df %>% filter(FDR < 0.05, est < 0) %>% arrange(est)        %>% slice_head(n = 3)
  bind_rows(up, down)
}
lab_ad <- pick_top3_labels(ad_df)
lab_k  <- pick_top3_labels(k_df)

# AD volcano 
ad_vol <- ggplot(ad_df, aes(x = est, y = mlog10FDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = FDR_LINE_SIZE, color = FDR_LINE_COLOR) +
  geom_vline(xintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_point(aes(color = sig_dir, alpha = alpha_val), size = 1.8, stroke = 0) +
  scale_color_manual(values = pal_ad, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  scale_y_continuous(limits = SIG_LIMS, breaks = SIG_BRKS, expand = expansion(mult = 0)) +
  ggrepel::geom_text_repel(
    data = lab_ad, aes(label = clean_label, color = sig_dir),
    size = 3, box.padding = 0.35, point.padding = 0.25, max.overlaps = 1000,
    force = 1, segment.alpha = 0.6, segment.size = 0.3, seed = 123, show.legend = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    plot.title      = element_blank(),
    plot.margin     = m_tight
  )

# K3–K1 volcano (rotated left) 
k_vol_rot <- ggplot(k_df, aes(x = est, y = mlog10FDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             linewidth = FDR_LINE_SIZE, color = FDR_LINE_COLOR) +
  geom_vline(xintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_point(aes(color = sig_dir, alpha = alpha_val), size = 1.8, stroke = 0) +
  scale_color_manual(values = pal_k, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  coord_flip() +
  scale_y_reverse(limits = rev(SIG_LIMS), breaks = rev(SIG_BRKS), expand = expansion(mult = 0)) +
  ggrepel::geom_text_repel(
    data = lab_k, aes(label = clean_label, color = sig_dir),
    size = 3, box.padding = 0.35, point.padding = 0.25, max.overlaps = 1000,
    force = 1, segment.alpha = 0.6, segment.size = 0.3, seed = 123, show.legend = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    plot.title      = element_blank(),
    plot.margin     = m_tight
  )

# Scatter
sc_df <- scatter_comp$data %>%
  mutate(alpha_val = ifelse(sig_cat == "ns", ALPHA_NS, ALPHA_SIG))

# labels
label_idx <- rep(FALSE, nrow(sc_df))
idx_sig   <- which(sc_df$sig_cat != "ns")
ord_sig   <- idx_sig[order(-sc_df$abs_sum[idx_sig])]
keep <- integer(0)
for (q in c("RU","LU","RD","LD")) {
  qidx <- ord_sig[sc_df$quadrant[ord_sig] == q]
  if (length(qidx)) keep <- c(keep, head(qidx, 2))
}
remain <- setdiff(ord_sig, keep)
keep <- unique(c(keep, head(remain[order(-sc_df$abs_sum[remain])], max(0, 25 - length(keep)))))
idx_both <- which(sc_df$sig_cat == "both")
if (length(idx_both) >= 1) keep <- unique(c(keep, idx_both[which.max(sc_df$abs_sum[idx_both])]))
label_idx[keep] <- TRUE

scatter_p <- ggplot(sc_df, aes(x = est_x, y = est_y, color = sig_cat, alpha = alpha_val)) +
  geom_hline(yintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_vline(xintercept = 0, linetype = "solid",
             linewidth = ZERO_LINE_SIZE, color = ZERO_LINE_COLOR) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = DIAG_LINE_SIZE, color = DIAG_LINE_COLOR) +
  geom_point(size = 1.7, stroke = 0) +
  scale_color_manual(values = pal_sc, breaks = names(pal_sc), drop = FALSE, guide = "none") +
  scale_alpha_identity() +
  scale_x_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-LIM_FC, LIM_FC), breaks = BRKS_FC, expand = expansion(mult = 0)) +
  coord_equal() +
  ggrepel::geom_text_repel(
    data = sc_df[label_idx, , drop = FALSE],
    aes(label = label_plot),
    size = 3, box.padding = 0.4, point.padding = 0.3,
    force = 2.0, min.segment.length = 0, segment.alpha = 0.6, segment.size = 0.3, seed = 123
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title      = element_blank(),
    plot.title      = element_blank(),
    plot.margin     = m_tight
  )

# Layout
g_combo <- ggdraw() +
  draw_plot(k_vol_rot, x = 0/3, y = 0/3, width = 1/3, height = 2/3) +
  draw_plot(ad_vol,    x = 1/3, y = 2/3, width = 2/3, height = 1/3) +
  draw_plot(scatter_p, x = 1/3, y = 0/3, width = 2/3, height = 2/3)

# Save
out_pdf <- file.path(OUT_DIR, "AD_K3K1_scatter_cowplot_grid_ticks_top3.pdf")
ggsave(out_pdf, plot = g_combo, device = "pdf", width = 7.5, height = 7.5, units = "in")
message("✅ Saved:", normalizePath(out_pdf))

























## ===== Residualize metabolites by covariates =================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
})

## ---------- Setups ----
# metabolites: data.frame with "projid" + metabolite columns, created above

COVARS <- c("pmi","study","fixation_interval","msex")

suppressPackageStartupMessages({
  library(dplyr); library(rlang)
})
metabolites <- metabolites %>% mutate(projid = as.character(projid))
meta        <- meta %>% mutate(projid = as.character(projid))

## ---------- 1) Build covariate df ----

stopifnot(all(COVARS %in% names(meta)))

cov_df <- meta %>%
  select(projid, all_of(COVARS)) %>%
  distinct(projid, .keep_all = TRUE) %>%         
  mutate(
    # numeric covariates
    pmi               = suppressWarnings(as.numeric(pmi)),
    fixation_interval = suppressWarnings(as.numeric(fixation_interval)),
    # categorical covariates
    study             = as.factor(study),
    msex              = as.factor(msex)
  )

## ---------- 2) Join covariates to metabolites; keep only rows with complete covariates ----
dat <- metabolites %>%
  left_join(cov_df, by = "projid") %>%
  filter(if_all(all_of(COVARS), ~ !is.na(.)))     

# Identify metabolite columns 
meta_cols <- setdiff(names(metabolites), "projid")
stopifnot(length(meta_cols) > 0)

## ---------- 3) Residualization function (y ~ covars per metabolite) ----
residualize_one <- function(y, X_df) {
  # rows with complete covariates and non-missing metabolite
  ok <- complete.cases(X_df) & !is.na(y)
  res <- rep(NA_real_, length(y))
  coefs <- NA
  
  if (sum(ok) >= (ncol(model.matrix(~ ., data = X_df)) + 1)) {
    # Build design matrix with intercept and default contrasts
    X <- model.matrix(~ ., data = X_df[ok, , drop = FALSE])
    fit <- lm.fit(x = X, y = y[ok])
    res[ok] <- y[ok] - (X %*% fit$coefficients)[, 1]
    coefs <- setNames(unname(fit$coefficients), colnames(X))
  }
  list(residuals = res, coefs = coefs)
}

# Prepare the covariate design table (no ID)
X_all <- dat %>% select(all_of(COVARS))

## ---------- 4) Apply to every metabolite ----
res_list <- lapply(meta_cols, function(met) {
  out <- residualize_one(dat[[met]], X_all)
  list(name = met, resid = out$residuals, coefs = out$coefs)
})

# Collect residuals into a wide data.frame
res_mat <- do.call(cbind, lapply(res_list, `[[`, "resid"))
colnames(res_mat) <- paste0(meta_cols, "_res")
metabolites_resid <- cbind(projid = dat$projid, as.data.frame(res_mat))

## ---------- 5) Save outputs ----
dir.create("Results/Metabolites_residualized", recursive = TRUE, showWarnings = FALSE)
saveRDS(metabolites_resid, file = "Results/Metabolites_residualized/metabolites_residualized_by_pmi_study_fix_msex.rds")

coef_tbl <- purrr::map_dfr(res_list, function(x) {
  tibble(metabolite = x$name,
         term = names(x$coefs),
         estimate = unname(x$coefs))
})
saveRDS(coef_tbl, file = "Results/Metabolites_residualized/metabolite_models_coefficients.rds")

## ---------- 6) sanity checks ----
message("Residualized matrix dims: ", nrow(metabolites_resid), " x ", ncol(metabolites_resid))
message("Example head (first 5 residualized metabolites):")
print(metabolites_resid[1:5, c("projid", colnames(metabolites_resid)[2:min(6, ncol(metabolites_resid))])])






## ===== Projid overlap metabolites and single cell data

meta<-load_BIG_meta_data()

library(eulerr)

ids <- list(
  metabolites_resid = unique(as.character(metabolites_resid$projid)),
  merged_df_final   = unique(as.character(merged_df_final$projid)),
  meta              = unique(as.character(meta$projid))
)

plot(
  euler(ids),
  quantities = TRUE, labels = TRUE,
  fills = list(fill = c("#4E79A7","#F28E2B","#E15759"), alpha = 0.35),
  edges = list(col = "grey20"),
  main = "projid overlap"
)

# Define ID sets
ids <- list(
  metabolites_resid = unique(as.character(metabolites_resid$projid)),
  merged_df_final   = unique(as.character(merged_df_final$projid))
)

# Generate Euler diagram
pdf("Results/projid_overlap_metabolites_vs_merged_df_final.pdf", width = 4, height = 4)

plot(
  euler(ids),
  quantities = TRUE,
  labels = TRUE,
  fills = list(fill = c("#606c38", "#bc6c25"), alpha = 0.85),
  edges = list(col = "grey20"),
  main = "projid overlap (metabolites vs merged_df_final)"
)

dev.off()

# --- Collect unique IDs ---
met_ids <- unique(as.character(metabolites_resid$projid))
sc_ids  <- unique(as.character(merged_df_final$projid))

# --- Compute sets ---
met_only      <- setdiff(met_ids, sc_ids)
singlecell_only <- setdiff(sc_ids, met_ids)
both          <- intersect(met_ids, sc_ids)

# --- Assemble into a tidy data frame ---
projid_df <- data.frame(
  projid = c(met_only, singlecell_only, both),
  group  = c(
    rep("metabolites_only", length(met_only)),
    rep("singlecell_only", length(singlecell_only)),
    rep("metabolites_and_singlecell", length(both))
  ),
  stringsAsFactors = FALSE
)

head(projid_df)


## ===== final GLMNET approach =====
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2)
  library(pROC); library(yardstick); library(glmnet); library(broom)
  library(patchwork)
})

set.seed(123)
OUT <- "Results/GLMNET_AD_solid"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---------- 0) Inputs & labels from meta  +  collapse replicates by projid (median) ----------
suppressPackageStartupMessages({ library(dplyr) })

meta <- load_BIG_meta_data()

stopifnot(all(c("projid") %in% names(metabolites_resid)))
stopifnot(all(c("projid","AD_status") %in% names(meta)))
stopifnot(length(met_only) > 0, length(both) > 0)

# 1) Build a deduplicated label table (one row per projid)
label_raw <- meta %>%
  mutate(projid = as.character(projid)) %>%
  select(projid, AD_status) %>%
  filter(!is.na(AD_status))

# check for conflicts
label_chk <- label_raw %>%
  group_by(projid) %>%
  summarise(n = n(), n_pos = sum(AD_status == 1), n_neg = sum(AD_status == 0), .groups = "drop")

conflicts <- label_chk %>% filter(n_pos > 0 & n_neg > 0)
if (nrow(conflicts) > 0) {
  warning("Conflicting AD_status for some projid in meta; resolving by majority vote (ties -> 1).")
}

label_df <- label_raw %>%
  group_by(projid) %>%
  summarise(
    AD_status = {
      v <- AD_status
      ones <- sum(v == 1); zeros <- sum(v == 0)
      if (zeros > ones) 0L else 1L
    },
    .groups = "drop"
  ) %>%
  mutate(AD_status = as.integer(AD_status))

# 2) Collapse metabolite replicates (median)
agg_fun <- function(x) suppressWarnings(stats::median(x, na.rm = TRUE))
feat_cols_all <- grep("_res$", names(metabolites_resid), value = TRUE)
stopifnot(length(feat_cols_all) > 0)

met_agg <- metabolites_resid %>%
  mutate(projid = as.character(projid)) %>%
  group_by(projid) %>%
  summarise(
    dplyr::across(all_of(feat_cols_all), agg_fun),
    n_reps = dplyr::n(),
    .groups = "drop"
  )

# 3) Join (labels are now unique per projid)
dat_all <- met_agg %>%
  inner_join(label_df, by = "projid") %>%
  filter(AD_status %in% c(0,1))

feat_cols <- grep("_res$", names(dat_all), value = TRUE)
stopifnot(length(feat_cols) > 0)

# 4) Split & matrices
train_df <- dat_all %>% filter(projid %in% as.character(met_only))
test_df  <- dat_all %>% filter(projid %in% as.character(both))


stopifnot(nrow(train_df) == dplyr::n_distinct(train_df$projid))
stopifnot(nrow(test_df)  == dplyr::n_distinct(test_df$projid))

X_tr_raw <- as.matrix(train_df[, feat_cols, drop = FALSE])
X_te_raw <- as.matrix(test_df[,  feat_cols, drop = FALSE])
y_tr <- train_df$AD_status
y_te <- test_df$AD_status

message("Projids in train_df:", length(train_df$projid), " Projids in test_df:", length(test_df$projid))

message("Replicate collapse complete: TRAIN projids = ", nrow(train_df),
        " | TEST projids = ", nrow(test_df))



## ---------- 1) Train-only imputation & scaling  --------
med_tr <- apply(X_tr_raw, 2, function(v) suppressWarnings(median(v, na.rm = TRUE)))
med_tr[!is.finite(med_tr)] <- 0
impute <- function(M, med) { if (anyNA(M)) { idx <- which(is.na(M), arr.ind = TRUE); M[idx] <- med[idx[,2]] }; M }
X_tr_imp <- impute(X_tr_raw, med_tr); X_te_imp <- impute(X_te_raw, med_tr)

mu  <- colMeans(X_tr_imp)
sdv <- apply(X_tr_imp, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
scale_mat <- function(M, mu, sdv) sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
X_tr_sc <- scale_mat(X_tr_imp, mu, sdv); X_te_sc <- scale_mat(X_te_imp, mu, sdv)

## ---------- 2) Nested CV: tune alpha & lambda --------
# Outer CV selects alpha; inner CV selects lambda for each alpha.
alphas <- c(0.0, 0.25, 0.5, 0.75, 1.0)
Kouter <- 5
fold_id <- sample(rep(1:Kouter, length.out = length(y_tr)))

w_class <- if (sum(y_tr==1)>0) ifelse(y_tr==1, sum(y_tr==0)/sum(y_tr==1), 1) else rep(1, length(y_tr))

outer_auc <- data.frame(alpha = numeric(), auc = numeric())
for (a in alphas) {
  cv_inner <- cv.glmnet(X_tr_sc, y_tr, family = "binomial", alpha = a,
                        type.measure = "auc", weights = w_class, nfolds = 10,
                        standardize = FALSE)
  lam <- cv_inner$lambda.1se
  # outer CV estimate using chosen lam
  preds <- rep(NA_real_, length(y_tr))
  for (k in 1:Kouter) {
    tr_idx <- which(fold_id != k); va_idx <- which(fold_id == k)
    fit_k <- glmnet(X_tr_sc[tr_idx, , drop = FALSE], y_tr[tr_idx],
                    family = "binomial", alpha = a, lambda = lam,
                    weights = w_class[tr_idx], standardize = FALSE)
    preds[va_idx] <- as.numeric(predict(fit_k, X_tr_sc[va_idx, , drop = FALSE], type = "response"))
  }
  auc_k <- as.numeric(pROC::auc(y_tr, preds))
  outer_auc <- rbind(outer_auc, data.frame(alpha = a, auc = auc_k))
}
best_alpha <- outer_auc$alpha[which.max(outer_auc$auc)]
readr::write_csv(outer_auc, file.path(OUT, "nested_outer_auc_alpha.csv"))

## ---------- 3) Fit final model on all TRAIN with best alpha + lambda.1se (inner CV) ----
cv_final <- cv.glmnet(X_tr_sc, y_tr, family = "binomial", alpha = best_alpha,
                      type.measure = "auc", weights = w_class, nfolds = 10, standardize = FALSE)
lam_final <- cv_final$lambda.1se
fit_glm <- glmnet(X_tr_sc, y_tr, family = "binomial", alpha = best_alpha,
                  lambda = lam_final, weights = w_class, standardize = FALSE)

## ---------- 4) Predict on HELD-OUT TEST --------
p_te <- as.numeric(predict(fit_glm, X_te_sc, type = "response"))

## ---------- 5) Metrics, threshold grid, and plots --------
brier <- mean((p_te - y_te)^2)
logloss <- { p <- pmin(pmax(p_te,1e-15),1-1e-15); -mean(y_te*log(p)+(1-y_te)*log(1-p)) }
auc_te <- as.numeric(pROC::auc(y_te, p_te))
eval_te <- tibble(AD_status = factor(y_te, levels=c(0,1), labels=c("nonAD","AD")), .pred_AD = p_te)
pr_te <- yardstick::pr_auc(eval_te, truth = AD_status, .pred_AD)$.estimate
roc_te <- yardstick::roc_curve(eval_te, truth = AD_status, .pred_AD)
prc_te <- yardstick::pr_curve(eval_te, truth = AD_status, .pred_AD)

# Youden threshold
th_te <- roc_te %>%
  mutate(J = sensitivity + specificity - 1) %>%
  slice_max(J, n = 1) %>% slice_head(n = 1) %>% pull(.threshold) %>% as.numeric()
if (!is.finite(th_te)) th_te <- 0.5

# Threshold grid
ths <- seq(0,1,by=0.01)
grid <- lapply(ths, function(t) {
  lab <- factor(ifelse(p_te >= t, "AD","nonAD"), levels=c("nonAD","AD"))
  tibble(
    threshold = t,
    accuracy = yardstick::accuracy_vec(eval_te$AD_status, lab),
    sensitivity = yardstick::sensitivity_vec(eval_te$AD_status, lab),
    specificity = yardstick::specificity_vec(eval_te$AD_status, lab),
    ppv = yardstick::ppv_vec(eval_te$AD_status, lab),
    npv = yardstick::npv_vec(eval_te$AD_status, lab)
  )
}) %>% bind_rows()
write_csv(grid, file.path(OUT, "threshold_grid_test.csv"))

# Save metrics summary
write_csv(
  tibble(model="GLMNET", alpha=best_alpha, lambda=lam_final,
         auc=auc_te, pr_auc=pr_te, brier=brier, log_loss=logloss,
         threshold_Youden=th_te),
  file.path(OUT,"metrics_test_summary.csv")
)

# Plots
p_roc <- ggplot(roc_te, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() + coord_equal() + theme_classic() + ggtitle(sprintf("ROC (AUC = %.3f)", auc_te))
ggsave(file.path(OUT,"ROC_glmnet.pdf"), p_roc, width=5, height=4)

p_pr <- ggplot(prc_te, aes(x = recall, y = precision)) +
  geom_path() + theme_classic() + ggtitle(sprintf("PR (PR-AUC = %.3f)", pr_te))
ggsave(file.path(OUT,"PR_glmnet.pdf"), p_pr, width=5, height=4)


p_dens <- ggplot(
  tibble(pred = p_te, AD_status = factor(y_te, levels = c(0, 1), labels = c("nonAD", "AD"))),
  aes(x = pred, fill = AD_status)
) +
  geom_density(alpha = 0.45) +
  geom_vline(xintercept = th_te, linetype = "dashed") +
  scale_fill_manual(values = c("nonAD" = "blue", "AD" = "red")) +
  theme_classic() +
  labs(
    x = "Predicted AD probability",
    y = "Density",
    title = "Prediction distribution"
  )
ggsave(file.path(OUT, "pred_density_by_true_status.pdf"), p_dens, width = 6, height = 4)

# Calibration (reliability) curve via deciles
cal_df <- tibble(pred=p_te, y=y_te) %>%
  mutate(bin = cut(pred, breaks = quantile(pred, probs = seq(0,1,0.1), na.rm=TRUE),
                   include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(pred_mean = mean(pred), obs_rate = mean(y), .groups="drop")
p_cal <- ggplot(cal_df, aes(x=pred_mean, y=obs_rate)) +
  geom_abline(slope=1, intercept=0, linetype="dotted") +
  geom_point() + geom_line() + coord_equal() +
  theme_classic() + labs(x="Mean predicted probability", y="Observed AD rate", title="Calibration (deciles)")
ggsave(file.path(OUT,"calibration_deciles.pdf"), p_cal, width=5.2, height=5)

## ---------- 6) Export per-donor TEST predictions --------
pred_test <- tibble(
  projid = test_df$projid,
  AD_status_true = factor(y_te, levels=c(0,1), labels=c("nonAD","AD")),
  pred_AD_prob_glmnet = p_te,
  pred_AD_label_glmnet = factor(ifelse(p_te >= th_te, "AD","nonAD"), levels=c("nonAD","AD")),
  threshold_used = th_te
)
write_csv(pred_test, file.path(OUT,"predictions_test_both.csv"))

## ---------- 7) merge into merged_df_final
if (exists("merged_df_final")) {
  merged_df_final <- merged_df_final %>%
    mutate(projid = as.character(projid)) %>%
    left_join(pred_test %>% select(projid, pred_AD_prob_glmnet, pred_AD_label_glmnet),
              by = "projid")
  saveRDS(merged_df_final, file.path(OUT,"merged_df_final_with_glmnet_preds.rds"))
}

## ---------- 8) Feature importance --------
# (A) Coefficients at final lambda (standardized scale)
coef_tbl <- as.matrix(coef(fit_glm))[,1, drop = FALSE]
coef_tbl <- tibble(Feature = rownames(coef_tbl), Coef = as.numeric(coef_tbl)) %>%
  filter(Feature != "(Intercept)") %>% arrange(desc(abs(Coef)))
write_csv(coef_tbl, file.path(OUT,"coefficients_final.csv"))

p_coef <- ggplot(head(coef_tbl, 25), aes(x = reorder(Feature, abs(Coef)), y = Coef)) +
  geom_col() + coord_flip() + theme_classic() + labs(x=NULL, y="Standardized coefficient",
                                                     title="Top 25 |β| (GLMNET)")
ggsave(file.path(OUT,"importance_coefficients_top25.pdf"), p_coef, width=6.5, height=5)

# (B) Permutation importance on TEST (AUC drop per feature)
perm_imp <- lapply(feat_cols, function(f) {
  X_perm <- X_te_sc
  X_perm[, f] <- sample(X_perm[, f])
  p_perm <- as.numeric(predict(fit_glm, X_perm, type="response"))
  data.frame(Feature=f, delta_auc = auc_te - as.numeric(pROC::auc(y_te, p_perm)))
}) %>% bind_rows() %>% arrange(desc(delta_auc))
write_csv(perm_imp, file.path(OUT,"importance_permutation_auc_drop_test.csv"))

p_perm <- ggplot(head(perm_imp, 25), aes(x = reorder(Feature, delta_auc), y = delta_auc)) +
  geom_col() + coord_flip() + theme_classic() + labs(x=NULL, y="AUC drop on permutation",
                                                     title="Permutation importance (TEST)")
ggsave(file.path(OUT,"importance_permutation_top25.pdf"), p_perm, width=6.5, height=5)

# (C) Stability selection (bootstrap selection frequency on TRAIN)
B <- 100
set.seed(123)
sel_freq <- setNames(numeric(length(feat_cols)), feat_cols)
for (b in 1:B) {
  idx <- sample.int(nrow(X_tr_sc), replace = TRUE)
  cvb <- cv.glmnet(X_tr_sc[idx, , drop = FALSE], y_tr[idx],
                   family="binomial", alpha=best_alpha, type.measure="auc",
                   weights=w_class[idx], nfolds=5, standardize=FALSE)
  fitb <- glmnet(X_tr_sc[idx, , drop = FALSE], y_tr[idx], family="binomial",
                 alpha=best_alpha, lambda=cvb$lambda.1se, weights=w_class[idx], standardize=FALSE)
  nz <- rownames(coef(fitb))[which(as.numeric(coef(fitb)) != 0)]
  nz <- setdiff(nz, "(Intercept)")
  sel_freq[nz] <- sel_freq[nz] + 1
}
stab_tbl <- tibble(Feature = names(sel_freq), sel_freq = as.numeric(sel_freq)/B) %>%
  arrange(desc(sel_freq))
write_csv(stab_tbl, file.path(OUT,"importance_stability_selection_freq.csv"))

p_stab <- ggplot(head(stab_tbl, 25), aes(x = reorder(Feature, sel_freq), y = sel_freq)) +
  geom_col() + coord_flip() + theme_classic() +
  labs(x=NULL, y="Selection frequency", title="Stability selection (TRAIN bootstraps)")
ggsave(file.path(OUT,"importance_stability_top25.pdf"), p_stab, width=6.5, height=5)

message("✅ GLMNET pipeline done. Outputs in: ", normalizePath(OUT))


