## ===== Start

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ---------- Lookup for DEG and pseudobulk files ----------

deg_name_lut <- c(
  "Mic.P2RY12"               = "Mic P2RY12",
  "Exc.L6.THEMIS.NFIA"       = "Exc L6 THEMIS NFIA",
  "Exc.L4.5.RORB.GABRG1"     = "Exc L4-5 RORB GABRG1",
  "Exc.L4.5.RORB.IL1RAPL2"   = "Exc L4-5 RORB IL1RAPL2",
  "Inh.LAMP5.NRG1..Rosehip." = "Inh LAMP5 NRG1 (Rosehip)",
  "Ast.GRM3"                 = "Ast GRM3",
  "Ast.DPP10"                = "Ast DPP10",
  "Oli"                      = "Oli",
  "OPCs"                     = "OPC"
)

pb_path_from_celltype <- function(ct_label, pb_dir) {
  file.path(pb_dir, paste0(ct_label, ".rds"))
}

find_deg_path <- function(ct_label, deg_dir) {
  if (!ct_label %in% names(deg_name_lut)) return(NA_character_)
  stem <- deg_name_lut[[ct_label]]
  path <- file.path(deg_dir, paste0(stem, ".rds"))
  if (file.exists(path)) path else NA_character_
}

## ---------------- palette & transforms ----------------

.z_cap <- function(M, cap = 2) {
  Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0
  if (!is.null(cap)) { Z[Z > cap] <- cap; Z[Z < -cap] <- -cap }
  Z
}

.make_palette <- function(low_col="#fefae0", mid_col=NULL, high_col="#606c38", n=201, cap_z=2) {
  brks <- if (is.null(cap_z)) seq(-3, 3, length.out = n) else seq(-cap_z, cap_z, length.out = n)
  cols <- if (is.null(mid_col)) colorRampPalette(c(low_col, high_col))(length(brks)-1)
  else                  colorRampPalette(c(low_col, mid_col, high_col))(length(brks)-1)
  list(cols = cols, breaks = brks)
}

.pad_to_order <- function(Z_present, donors_order, donors_present, na_val = NA_real_) {
  Z_full <- matrix(na_val, nrow = nrow(Z_present), ncol = length(donors_order),
                   dimnames = list(rownames(Z_present), donors_order))
  if (length(donors_present)) Z_full[, donors_present] <- Z_present
  Z_full
}

## ---------------- selection: top up/down from DESeqResults ----------------

select_top_genes_from_DESeq <- function(
    deg_rds_path,
    pb_gene_ids,
    deg_slot      = "subgroups",
    top_up        = 25,
    top_down      = 25,
    fdr_cut       = 0.05,
    lfc_cut       = 0.25,
    verbose       = TRUE
){
  obj <- readRDS(deg_rds_path)
  stopifnot(deg_slot %in% names(obj))
  res <- obj[[deg_slot]]
  df  <- as.data.frame(res)
  
  df$symbol <- if ("geneName" %in% names(df)) df$geneName else rownames(df)
  df <- df %>% filter(is.finite(padj), is.finite(log2FoldChange), symbol %in% pb_gene_ids)
  
  up   <- df %>% filter(padj < fdr_cut, log2FoldChange >=  lfc_cut) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n = top_up)   %>% mutate(dir = "UP")
  down <- df %>% filter(padj < fdr_cut, log2FoldChange <= -lfc_cut) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n = top_down) %>% mutate(dir = "DOWN")
  
  sel <- bind_rows(up, down)
  if (isTRUE(verbose)) message(sprintf("• Selected %d UP and %d DOWN (total %d).", nrow(up), nrow(down), nrow(sel)))
  sel %>% transmute(gene = symbol, dir, padj, log2FC = log2FoldChange)
}

## ---------------- plotting --------

plot_celltype_heatmap <- function(
    pb,                      # genes x donors (logCPM)
    sel_genes_tbl,
    donors_order,            # fixed order
    merged_df_final,         # for subgroup annotation
    subgroup_col   = "subgroup",
    cap_z          = 2,
    low_col        = "#fefae0",
    mid_col        = NULL,
    high_col       = "#606c38",
    na_col         = "white",
    row_kmeans_k   = NULL,   
    row_kmeans_nstart = 50,
    title          = NULL,
    show_rownames  = TRUE
){
  stopifnot(is.matrix(pb))
  donors_order <- as.character(donors_order)
  
  donors_present <- intersect(donors_order, colnames(pb))
  donors_missing <- setdiff(donors_order, donors_present)
  if (!length(donors_present)) stop("No overlap between donors_order and pseudobulk columns.")
  
  genes <- intersect(unique(sel_genes_tbl$gene), rownames(pb))
  if (length(genes) < 2) stop("Fewer than 2 selected genes exist in pseudobulk.")
  X <- pb[genes, donors_present, drop = FALSE]
  
  Z <- .z_cap(X, cap = cap_z)
  
  gaps_row <- NULL
  ann_row  <- NULL
  if (!is.null(row_kmeans_k) && row_kmeans_k > 1L && nrow(Z) >= row_kmeans_k) {
    set.seed(1)
    km <- kmeans(Z, centers = row_kmeans_k, nstart = row_kmeans_nstart)
    cl <- factor(km$cluster, levels = sort(unique(km$cluster)))
    cents <- sapply(levels(cl), function(k) colMeans(Z[cl == k, , drop = FALSE]))
    if (!is.matrix(cents)) cents <- matrix(cents, ncol = 1)
    cm <- suppressWarnings(cor(cents, use = "pairwise.complete.obs"))
    cm[!is.finite(cm)] <- 0; diag(cm) <- 1
    ord_lv <- hclust(as.dist(1 - cm), method = "ward.D2")$order
    cl <- factor(cl, levels = levels(cl)[ord_lv])
    .ord_rows <- function(M) {
      if (nrow(M) < 3) return(seq_len(nrow(M)))
      cm <- suppressWarnings(cor(t(M), use = "pairwise.complete.obs"))
      cm[!is.finite(cm)] <- 0; diag(cm) <- 1
      hclust(as.dist(1 - cm), method = "ward.D2")$order
    }
    ord <- integer(0)
    for (lv in levels(cl)) {
      idx <- which(cl == lv)
      ord <- c(ord, idx[.ord_rows(Z[idx, , drop = FALSE])])
    }
    Z  <- Z[ord, , drop = FALSE]
    cl <- cl[ord]
    ann_row  <- data.frame(Kmeans = cl, row.names = rownames(Z), check.names = FALSE)
    gaps_row <- cumsum(table(cl)); gaps_row <- gaps_row[-length(gaps_row)]
  }
  
  Z_full <- .pad_to_order(Z, donors_order, donors_present, na_val = NA_real_)
  
  ann_df <- merged_df_final %>%
    mutate(projid = as.character(projid)) %>%
    select(projid, !!rlang::sym(subgroup_col)) %>%
    distinct(projid, .keep_all = TRUE)
  ann_df <- ann_df[match(donors_order, ann_df$projid), , drop = FALSE]
  ann_col <- data.frame(subgroup = factor(ann_df[[subgroup_col]]),
                        row.names = donors_order, check.names = FALSE)
  
  pal <- .make_palette(low_col, mid_col, high_col, n = 201, cap_z = cap_z)
  if (is.null(title)) {
    title <- sprintf("Top DEGs (UP/DOWN) • donors fixed (%d present, %d missing)",
                     length(donors_present), length(donors_missing))
  }
  
  pheatmap(
    Z_full,
    color = pal$cols, breaks = pal$breaks, na_col = na_col,
    annotation_col = ann_col,
    annotation_row = ann_row,
    cluster_rows = FALSE, cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = show_rownames && (nrow(Z_full) <= 150),
    fontsize_row = 7,
    gaps_row = gaps_row,
    main = title
  )
  
  invisible(list(Z_present = Z, Z_full = Z_full,
                 donors_present = donors_present, donors_missing = donors_missing))
}

## ---------------- existence check  -----------

check_celltype_files <- function(cell_types,
                                 pseudobulk_dir = "Data/pseudobulk",
                                 deg_dir        = "Data/subgroup_DEGs",
                                 verbose        = TRUE) {
  df <- lapply(cell_types, function(ct) {
    pb_path  <- pb_path_from_celltype(ct, pseudobulk_dir)
    deg_path <- find_deg_path(ct, deg_dir)
    data.frame(
      cell_type       = ct,
      pseudobulk_ok   = file.exists(pb_path),
      deg_ok          = if (!is.na(deg_path)) file.exists(deg_path) else FALSE,
      pseudobulk_path = pb_path,
      deg_path        = deg_path %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  summary_df <- bind_rows(df)
  
  if (isTRUE(verbose)) {
    cat("\n=== Cell type file existence check ===\n")
    for (i in seq_len(nrow(summary_df))) {
      ct <- summary_df$cell_type[i]
      pb_ok  <- if (summary_df$pseudobulk_ok[i]) "✅" else "❌"
      deg_ok <- if (summary_df$deg_ok[i])        "✅" else "❌"
      cat(sprintf("%-24s | PB: %s  DEG: %s\n", ct, pb_ok, deg_ok))
    }
    cat("--------------------------------------\n")
    cat(sprintf("Total pseudobulk found: %d/%d\n",
                sum(summary_df$pseudobulk_ok), nrow(summary_df)))
    cat(sprintf("Total DEG found:        %d/%d\n",
                sum(summary_df$deg_ok), nrow(summary_df)))
  }
  invisible(summary_df)
}

## ---------------- runner -----------
run_for_celltypes <- function(
    cell_types,
    merged_df_final,
    donors_order,                       
    deg_slot         = "subgroups",
    top_up           = 20,
    top_down         = 20,
    fdr_cut          = 0.05,
    lfc_cut          = 0.25,
    cap_z            = 2,
    low_col          = "#606C38",
    mid_col          = "#dad7cd",
    high_col         = "#BC6C25",
    row_kmeans_k     = 2,
    pdf_dir          = NULL,
    pseudobulk_dir   = "Data/pseudobulk",
    deg_dir          = "Data/subgroup_DEGs"
){
  status <- check_celltype_files(cell_types, pseudobulk_dir, deg_dir, verbose = FALSE)
  keep_idx <- which(status$pseudobulk_ok & status$deg_ok)
  if (!length(keep_idx)) stop("No cell types with BOTH pseudobulk and DEG files. Aborting.")
  if (length(keep_idx) < length(cell_types)) {
    missing <- setdiff(cell_types, status$cell_type[keep_idx])
    message("⚠ Skipping (missing one or both files): ", paste(missing, collapse = ", "))
  }
  
  if (!is.null(pdf_dir)) dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
  
  out <- vector("list", length(keep_idx)); names(out) <- status$cell_type[keep_idx]
  for (ii in seq_along(keep_idx)) {
    i  <- keep_idx[ii]
    ct <- status$cell_type[i]
    pb_path  <- status$pseudobulk_path[i]
    deg_path <- status$deg_path[i]
    
    message("=== ", ct, " ===")
    pb <- readRDS(pb_path)                   
    stopifnot(is.matrix(pb))
    
    sel <- select_top_genes_from_DESeq(
      deg_rds_path = deg_path,
      pb_gene_ids  = rownames(pb),
      deg_slot     = deg_slot,
      top_up       = top_up,
      top_down     = top_down,
      fdr_cut      = fdr_cut,
      lfc_cut      = lfc_cut,
      verbose      = TRUE
    )
    
    ttl <- sprintf("%s • top %d up + %d down (FDR<%.2g, |LFC|≥%.2g)",
                   ct, top_up, top_down, fdr_cut, lfc_cut)
    
    if (is.null(pdf_dir)) {
      plot_celltype_heatmap(
        pb              = pb,
        sel_genes_tbl   = sel,
        donors_order    = donors_order,
        merged_df_final = merged_df_final,
        subgroup_col    = "subgroup",
        cap_z           = cap_z,
        low_col         = low_col,
        mid_col         = mid_col,
        high_col        = high_col,
        row_kmeans_k    = row_kmeans_k,
        title           = ttl
      )
    } else {
      pdf_file <- file.path(pdf_dir, paste0(ct, ".DEG_heatmap.pdf"))
      grDevices::pdf(pdf_file, width = 10, height = 3, useDingbats = FALSE)
      plot_celltype_heatmap(
        pb              = pb,
        sel_genes_tbl   = sel,
        donors_order    = donors_order,
        merged_df_final = merged_df_final,
        subgroup_col    = "subgroup",
        cap_z           = cap_z,
        low_col         = low_col,
        mid_col         = mid_col,
        high_col        = high_col,
        row_kmeans_k    = row_kmeans_k,
        title           = ttl
      )
      grDevices::dev.off()
      message("Saved: ", normalizePath(pdf_file))
    }
    
    out[[ct]] <- sel
  }
  invisible(out)
}

## ---------------- usage ---------------- 

donor_order <- res_all$col_order

cell_types <- c("Mic.P2RY12","Exc.L6.THEMIS.NFIA","Exc.L4.5.RORB.GABRG1","Exc.L4.5.RORB.IL1RAPL2",
                "Inh.LAMP5.NRG1..Rosehip.","Ast.GRM3","Ast.DPP10","Oli","OPCs")

summary_df <- check_celltype_files(
  cell_types     = cell_types,
  pseudobulk_dir = "Data/pseudobulk",
  deg_dir        = "Data/subgroup_DEGs",
  verbose        = TRUE
)

run_for_celltypes(
  cell_types       = cell_types,
  merged_df_final  = merged_df_final,
  donors_order     = donor_order,
  deg_slot         = "subgroups",
  top_up           = 10,
  top_down         = 10,
  fdr_cut          = 1e-4,
  lfc_cut          = 0.25,
  cap_z            = 1,
  low_col          = "#606C38", mid_col = "#dad7cd", high_col = "#BC6C25",
  row_kmeans_k     = 2,
  pdf_dir          = "Results/Celltype_DEG_heatmaps",
  pseudobulk_dir   = "Data/pseudobulk",
  deg_dir          = "Data/subgroup_DEGs"
)
