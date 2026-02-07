## ==== 1) Volcanoes for DESeqResults ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# ---- load DESeq results
MicP2RY12_DEGs<-readRDS('Data/subgroup_DEGs/Mic P2RY12.rds')

# ---- paths
OUT_BASE  <- "Results_subgroup_DEGs"
PLOT_DIR  <- file.path(OUT_BASE, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Volcano helpers 

# -- tidy 
deseq_to_df_volcano <- function(res) {
  df <- as.data.frame(res)
  gene <- if ("geneName" %in% names(df)) df$geneName else rownames(df)
  pval <- suppressWarnings(as.numeric(df$pvalue))
  padj <- suppressWarnings(as.numeric(df$padj))
  FDR  <- if ("pvalue" %in% names(df) && any(is.finite(pval))) {
    p.adjust(pval, method = "fdr")
  } else if ("padj" %in% names(df)) {
    padj
  } else {
    rep(NA_real_, nrow(df))
  }
  out <- data.frame(
    gene     = as.character(gene),
    log2FC   = suppressWarnings(as.numeric(df$log2FoldChange)),
    pval     = pval,
    padj     = padj,
    FDR      = as.numeric(FDR),
    stringsAsFactors = FALSE
  )
  out$mlog10FDR <- -log10(pmax(out$FDR, .Machine$double.xmin))
  out
}

make_deseq_volcano <- function(df,
                               title,
                               colors       = list(pos = "#D55E00", neg = "#0072B2", ns = "grey70"),
                               point_size   = 1.6,
                               point_alpha  = 0.85,
                               xlim_fixed   = c(-3.5, 3.5),
                               ylim_fixed   = c(0, 35),
                               fdr_thr      = 0.05,
                               lfc_thr      = 1.0,
                               label_per_side = 20,
                               fdr_fill2    = 0.20,
                               seed         = 123) {
  if (!"mlog10FDR" %in% names(df) && "FDR" %in% names(df)) {
    df$mlog10FDR <- -log10(pmax(df$FDR, .Machine$double.xmin))
  }
  
  d <- df
  is_sig <- is.finite(d$FDR) & d$FDR < fdr_thr
  d$sig_dir <- "ns"
  d$sig_dir[is_sig & d$log2FC > 0] <- "pos_sig"
  d$sig_dir[is_sig & d$log2FC < 0] <- "neg_sig"
  pal <- c(pos_sig = colors$pos, neg_sig = colors$neg, ns = colors$ns)
  
  p <- ggplot2::ggplot(d, aes(x = log2FC, y = mlog10FDR)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    ggplot2::geom_vline(xintercept = 0,               linetype = "solid",  linewidth = 0.5, color = "grey40") +
    ggplot2::geom_point(aes(color = sig_dir), alpha = point_alpha, size = point_size, stroke = 0) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::coord_cartesian(xlim = xlim_fixed, ylim = ylim_fixed, clip = "off") +
    ggplot2::labs(title = title, x = expression(log[2](FC)), y = expression(-log[10](FDR))) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0))
  
  pick_side <- function(side_df, n) {
    ord <- function(z) z[order(z$FDR, -abs(z$log2FC)), , drop = FALSE]
    t1 <- subset(side_df, FDR < fdr_thr & abs(log2FC) >= lfc_thr)
    t2 <- subset(side_df, FDR < fdr_thr)
    t3 <- subset(side_df, FDR < fdr_fill2)
    t4 <- side_df
    out <- ord(t1)
    if (nrow(out) < n) out <- unique(rbind(out, ord(t2)))
    if (nrow(out) < n) out <- unique(rbind(out, ord(t3)))
    if (nrow(out) < n) out <- unique(rbind(out, ord(t4)))
    head(out, n)
  }
  
  pool  <- d[!is.na(d$gene) & d$gene != "" & is.finite(d$FDR) & is.finite(d$log2FC), , drop = FALSE]
  up_df <- pick_side(subset(pool, log2FC > 0),  label_per_side)
  dn_df <- pick_side(subset(pool, log2FC < 0),  label_per_side)
  lab_df <- rbind(up_df, dn_df)
  lab_df$clean_gene <- gsub("\\.", " ", lab_df$gene)
  
  if (nrow(lab_df)) {
    set.seed(seed)
    p <- p + ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = clean_gene, color = sig_dir),
      size = 3, max.overlaps = 100, box.padding = 0.35, point.padding = 0.25,
      force = 0.9, min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6
    )
  }
  
  n_up   <- sum(d$sig_dir == "pos_sig", na.rm = TRUE)
  n_down <- sum(d$sig_dir == "neg_sig", na.rm = TRUE)
  message(sprintf("[volcano] %s — FDR<%.2f: up=%d, down=%d; labeled up=%d, down=%d",
                  title, fdr_thr, n_up, n_down, nrow(up_df), nrow(dn_df)))
  p
}

volcano_from_deseqres <- function(res,
                                  title,
                                  file_stub,
                                  colors    = list(pos = "#D55E00", neg = "#0072B2", ns = "grey70"),
                                  out_dir   = PLOT_DIR,
                                  ...) {
  df <- deseq_to_df_volcano(res)
  p  <- make_deseq_volcano(df, title = title, colors = colors, ...)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, ".png"))
  ggsave(pdf_file, plot = p, width = 6, height = 5, device = cairo_pdf)
  ggsave(png_file, plot = p, width = 6, height = 5, dpi = 300)
  message("[save] ", normalizePath(pdf_file))
  message("[save] ", normalizePath(png_file))
  invisible(p)
}



# ---- Build the three volcano plots

# 1) AD vs non-AD
volcano_from_deseqres(
  MicP2RY12_DEGs$casecontrol,
  title     = "Mic P2RY12: AD vs non-AD",
  file_stub = "Mic_P2RY12_AD_vs_nonAD",
  colors    = list(pos = "#F2696B", neg = "#869AD5", ns = "grey80")
)

# 2) K3 vs K1
volcano_from_deseqres(
  MicP2RY12_DEGs$cluster,
  title     = "Mic P2RY12: K3 vs K1",
  file_stub = "Mic_P2RY12_K3_vs_K1",
  colors    = list(pos = "#ad2524", neg = "#264653", ns = "grey80")
)

# 3) PathAD.TxAD vs PathNonAD.TxNonAD
volcano_from_deseqres(
  MicP2RY12_DEGs$subgroups,
  title     = "Mic P2RY12: PathAD.TxAD vs PathNonAD.TxNonAD",
  file_stub = "Mic_P2RY12_PathAD.TxAD_vs_PathNonAD.TxNonAD",
  colors    = list(pos = "#7B497D", neg = "#C8857F", ns = "grey75")
)



## ==== 2) TO ADD: scatter ----
## ==== 3) Rank-based pathway enrichment for DESeqResults (fgsea) ====
suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
})

OUT_BASE   <- "Results_subgroup_DEGs"
PATH_DIR   <- file.path(OUT_BASE, "pathways")
dir.create(PATH_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(msigdbr)
})

# Build a combined gene-set list for fgsea
build_gsets <- function(
    which = c("hallmark","kegg","reactome","go"),
    species = "Homo sapiens",
    prefix_source = TRUE
){
  # Map names to msigdbr keys
  coll_map <- list(
    hallmark = "H",
    kegg     = "C2:CP:KEGG",
    reactome = "C2:CP:REACTOME",
    go       = c("C5:BP","C5:CC","C5:MF")
  )
  want <- unlist(coll_map[intersect(which, names(coll_map))], use.names = FALSE)
  
  # Fetch per collection and bind
  dfl <- lapply(want, function(coll) {
    msigdbr(species = species, collection = coll) %>%
      dplyr::select(gs_name, gene_symbol, collection)
  })
  gs_tbl <- dplyr::bind_rows(dfl) %>% as.data.frame()
  
  # prefix set names with collection (helps when merging sources)
  set_names <- if (prefix_source && "collection" %in% names(gs_tbl)) {
    paste0(gs_tbl$collection, ":", gs_tbl$gs_name)
  } else {
    gs_tbl$gs_name
  }
  
  split(gs_tbl$gene_symbol, set_names)
}

# ==== Build gsets ====
gsets <- build_gsets(
  which = c("hallmark","kegg","reactome","go"),
  species = "Homo sapiens",
  prefix_source = TRUE
)

# sanity checks
length(gsets)
head(names(gsets), 8)


# ---- Helpers to tidy DESeqResults and build rank vector
deseq_to_df_fgsea <- function(res) {
  df <- as.data.frame(res)
  gene <- if ("geneName" %in% names(df)) df$geneName else rownames(df)
  pval <- suppressWarnings(as.numeric(df$pvalue))
  lfc  <- suppressWarnings(as.numeric(df$log2FoldChange))
  stat <- suppressWarnings(as.numeric(df$stat))
  FDR  <- if ("pvalue" %in% names(df) && any(is.finite(pval))) {
    p.adjust(pval, "fdr")
  } else if ("padj" %in% names(df)) as.numeric(df$padj) else rep(NA_real_, nrow(df))
  tibble::tibble(gene = as.character(gene),
                 log2FC = lfc,
                 pvalue = pval,
                 stat   = stat,
                 FDR    = FDR)
}

make_ranks <- function(df) {
  r <- df$stat
  if (!any(is.finite(r))) {
    r <- sign(df$log2FC) * -log10(pmax(df$pvalue, .Machine$double.xmin))
  }
  keep <- !is.na(df$gene) & df$gene != "" & is.finite(r)
  d <- df[keep, ]
  r <- r[keep]
  r <- tapply(r, d$gene, function(x) x[which.max(abs(x))])
  sort(unlist(r), decreasing = TRUE)
}

run_fgsea_on_res <- function(res, label,
                             gsets, nperm = 20000,
                             minSize = 15, maxSize = 500) {
  df    <- deseq_to_df_fgsea(res)
  ranks <- make_ranks(df)
  fg    <- fgsea::fgsea(pathways = gsets, stats = ranks,
                        nperm = nperm, minSize = minSize, maxSize = maxSize)
  fg <- dplyr::arrange(fg, padj) %>%
    dplyr::mutate(contrast = label,
                  direction = ifelse(NES >= 0, "up", "down")) %>%
    dplyr::select(contrast, pathway, NES, padj, size, leadingEdge, direction)
  out_csv <- file.path(PATH_DIR, paste0(gsub("[^A-Za-z0-9._-]", "_", label), "_fgsea.csv"))
  readr::write_csv(fg, out_csv)
  message("[fgsea] Saved: ", normalizePath(out_csv))
  fg
}

# ---- Run for three Mic P2RY12 contrasts
fg_case  <- run_fgsea_on_res(as.data.frame(MicP2RY12_DEGs$casecontrol), "AD_vs_nonAD", gsets)
fg_clust <- run_fgsea_on_res(as.data.frame(MicP2RY12_DEGs$cluster),     "K3_vs_K1",    gsets)
fg_sub   <- run_fgsea_on_res(as.data.frame(MicP2RY12_DEGs$subgroups),   "PathAD.TxAD_vs_PathNonAD.TxNonAD", gsets)

# ---- Merge NES across contrasts
all_fg <- bind_rows(fg_case, fg_clust, fg_sub)

# Keep only pathways significantly dysregulated in at least one contrast
fdr_thr <- 0.05
contrast_order <- c("AD_vs_nonAD","K3_vs_K1","PathAD.TxAD_vs_PathNonAD.TxNonAD")

sig_only <- all_fg %>%
  dplyr::filter(is.finite(padj) & padj < fdr_thr)

keep_paths <- unique(sig_only$pathway)

wide_nes <- all_fg %>%
  dplyr::filter(pathway %in% keep_paths) %>%
  dplyr::mutate(contrast = factor(contrast, levels = contrast_order)) %>%
  dplyr::select(pathway, contrast, NES, padj) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = c(NES, padj))

# Order by overall NES (desc)
nes_cols <- grep("^NES_", names(wide_nes), value = TRUE)
mean_abs_nes <- if (length(nes_cols)) rowMeans(abs(as.matrix(wide_nes[nes_cols])), na.rm = TRUE) else rep(NA_real_, nrow(wide_nes))
ord <- order(mean_abs_nes, decreasing = TRUE, na.last = TRUE)
wide_nes <- wide_nes[ord, , drop = FALSE]

message(sprintf("[wide_nes] kept %d pathways (FDR<%.2g in ≥1 contrast).", nrow(wide_nes), fdr_thr))

## ==== 4) Heatmap of significant pathway NES,non-significant entries white (padj>=0.05) ====
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
})

# matrices from wide_nes
nes_mat  <- as.matrix(wide_nes %>% select(starts_with("NES_")))
padj_mat <- as.matrix(wide_nes %>% select(starts_with("padj_")))
colnames(nes_mat)  <- sub("^NES_",  "", colnames(nes_mat))
colnames(padj_mat) <- sub("^padj_", "", colnames(padj_mat))
if ("pathway" %in% names(wide_nes)) rownames(nes_mat) <- wide_nes$pathway

# filter to rows significant in ≥1 contrast
fdr_cutoff <- 0.05
keep <- apply(padj_mat, 1, function(x) any(is.finite(x) & x < fdr_cutoff))
n_before <- nrow(padj_mat)
nes_mat  <- nes_mat[keep, , drop = FALSE]
padj_mat <- padj_mat[keep, , drop = FALSE]
message(sprintf("[filter] pathways kept: %d / %d (FDR < %.2f in ≥1 contrast)",
                nrow(nes_mat), n_before, fdr_cutoff))
stopifnot(nrow(nes_mat) > 0)

# color scale for NES
col_fun <- colorRamp2(c(-3, -1.5, 0, 1.5, 3),
                      c("#3957b7", "#9bb3e0", "white", "#f3a6a0", "#c43a31"))

# make non-significant cells white by setting them to NA
plot_nes <- nes_mat
plot_nes[padj_mat >= fdr_cutoff | is.na(padj_mat)] <- NA

# clustering
k_row_clusters <- 12
hc  <- hclust(dist(plot_nes %>% replace(is.na(.), 0)), method = "ward.D2") # replace NAs for dist
dnd <- as.dendrogram(hc)
cl_int <- cutree(hc, k = k_row_clusters)
names(cl_int) <- rownames(plot_nes)

# visual order 
row_ord_idx <- order.dendrogram(dnd)
rows_in_visual_order     <- rownames(plot_nes)[row_ord_idx]
clusters_in_visual_order <- cl_int[rows_in_visual_order]
ordered_cluster_ids      <- unique(clusters_in_visual_order)
k_eff <- length(ordered_cluster_ids)

# color 
max_brewer <- brewer.pal.info[brewer.pal.info$category == 'qual', 'maxcolors']
largest_set <- rownames(brewer.pal.info[brewer.pal.info$maxcolors == max(max_brewer), ])
base_colors <- brewer.pal(max(max_brewer), largest_set[1])
palette_hex <- if (k_eff > length(base_colors)) {
  colorRampPalette(base_colors)(k_eff)
} else base_colors[seq_len(k_eff)]

# visual names + legend labels
vis_names      <- paste0("vis_clus_", seq_len(k_eff))
legend_labels  <- sprintf("%s (ID: %s)", vis_names, ordered_cluster_ids)
id_to_visname  <- setNames(vis_names, ordered_cluster_ids)
id_to_label    <- setNames(legend_labels, ordered_cluster_ids)
label_to_color <- setNames(palette_hex, legend_labels)
id_to_color    <- setNames(palette_hex, ordered_cluster_ids)

# per-row cluster factor
cl_display_vec <- unname(id_to_label[as.character(cl_int)])
if (anyNA(cl_display_vec)) stop("[naming] Some cluster IDs not mapped to legend labels.")
cl_fac  <- factor(cl_display_vec, levels = legend_labels)
cl_cols <- setNames(palette_hex, levels(cl_fac))

# right-side annotation
row_anno <- rowAnnotation(
  Cluster = cl_fac,
  col = list(Cluster = cl_cols),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 9),
  width = unit(7, "mm")
)

# heatmap (white for NA cells)
ht <- Heatmap(
  plot_nes,
  name = "NES",
  col = col_fun,
  na_col = "white",                
  cluster_rows = dnd,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  heatmap_legend_param = list(direction = "horizontal"),
  row_dend_width = unit(3, "mm"),
  column_names_gp = gpar(fontsize = 9),
  border = TRUE
)

# draw + save
ht_list <- ht + row_anno
pdf(file.path(PATH_DIR, "Mic_P2RY12_pathway_NES_heatmap_visualClusters_whiteNA.pdf"),
    width = 1, height = 7.5, useDingbats = FALSE)
draw(ht_list, heatmap_legend_side = "bottom")
dev.off()

# CSV export
assignments <- tibble(
  pathway        = rownames(plot_nes),
  cluster_id     = unname(cl_int[rownames(plot_nes)]),
  visual_name    = unname(id_to_visname[as.character(cl_int[rownames(plot_nes)])]),
  legend_label   = unname(id_to_label[as.character(cl_int[rownames(plot_nes)])]),
  color_hex      = unname(id_to_color[as.character(cl_int[rownames(plot_nes)])]),
  visual_rank    = match(rownames(plot_nes), rows_in_visual_order)
)
csv_path <- file.path(PATH_DIR, "Mic_P2RY12_pathway_cluster_assignments_visualClusters_whiteNA.csv")
write_csv(assignments, csv_path)
message("[export] Cluster assignments: ", normalizePath(csv_path))

## ==== 5) Multi-cluster stacked barplots ====


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)   
})

# ---- constants: column & color mapping
.colmap <- list(
  casecontrol = list(padj="padj_AD_vs_nonAD",  nes="NES_AD_vs_nonAD",
                     pos_col="#F2696B", neg_col="#869AD5",
                     pos_lab="up in AD", neg_lab="up in non-AD"),
  cluster     = list(padj="padj_K3_vs_K1",     nes="NES_K3_vs_K1",
                     pos_col="#ad2524", neg_col="#264653",
                     pos_lab="up in K3", neg_lab="up in K1"),
  subgroups   = list(padj="padj_PathAD.TxAD_vs_PathNonAD.TxNonAD",
                     nes ="NES_PathAD.TxAD_vs_PathNonAD.TxNonAD",
                     pos_col="#7B497D", neg_col="#C8857F",
                     pos_lab="up in PathAD.TxAD", neg_lab="up in PathNonAD.TxNonAD")
)

.fill_levels <- c("casecontrol_pos","casecontrol_neg",
                  "cluster_pos","cluster_neg",
                  "subgroups_pos","subgroups_neg")

.fill_vals <- c(
  casecontrol_pos = .colmap$casecontrol$pos_col,
  casecontrol_neg = .colmap$casecontrol$neg_col,
  cluster_pos     = .colmap$cluster$pos_col,
  cluster_neg     = .colmap$cluster$neg_col,
  subgroups_pos   = .colmap$subgroups$pos_col,
  subgroups_neg   = .colmap$subgroups$neg_col
)

.fill_labs <- c(
  casecontrol_pos = .colmap$casecontrol$pos_lab,
  casecontrol_neg = .colmap$casecontrol$neg_lab,
  cluster_pos     = .colmap$cluster$pos_lab,
  cluster_neg     = .colmap$cluster$neg_lab,
  subgroups_pos   = .colmap$subgroups$pos_lab,
  subgroups_neg   = .colmap$subgroups$neg_lab
)

# ---- helper
.make_cluster_panel <- function(wide_nes, assignments,
                                cluster,            # e.g. "vis_clus_3" or numeric id or full legend
                                sort_by = c("casecontrol","cluster","subgroups"),
                                top_n   = 20,
                                xlim_max = NULL) {
  sort_by <- match.arg(sort_by)
  
  needed <- unlist(lapply(.colmap, `[`, c("padj","nes")))
  stopifnot(all(c("pathway", needed) %in% names(wide_nes)))
  
  # pick cluster rows
  cl_df <- assignments
  if (is.numeric(cluster)) {
    cl_df <- dplyr::filter(cl_df, cluster_id == !!cluster)
  } else if (grepl("^vis_clus_\\d+$", cluster)) {
    cl_df <- dplyr::filter(cl_df, visual_name == !!cluster)
  } else {
    cl_df <- dplyr::filter(cl_df, legend_label == !!cluster)
  }
  if (!nrow(cl_df)) stop("No pathways for cluster: ", cluster)
  title_lab <- unique(cl_df$legend_label)[1]
  
  # long format per pathway × contrast
  long_df <- dplyr::inner_join(
    cl_df[, c("pathway","cluster_id","visual_name","legend_label")],
    wide_nes[, c("pathway", needed)],
    by = "pathway"
  ) |>
    tidyr::pivot_longer(
      cols = -c(pathway, cluster_id, visual_name, legend_label),
      names_to = c(".value","contrast_key"),
      names_pattern = "^(padj|NES)_(.*)$"
    ) |>
    dplyr::mutate(
      contrast = dplyr::recode(contrast_key,
                               "AD_vs_nonAD" = "casecontrol",
                               "K3_vs_K1"    = "cluster",
                               "PathAD.TxAD_vs_PathNonAD.TxNonAD" = "subgroups"),
      direction = ifelse(NES >= 0, "pos", "neg"),
      fill_key  = paste(contrast, direction, sep = "_"),
      neglog10padj = -log10(pmax(padj, .Machine$double.xmin))
    ) |>
    dplyr::filter(contrast %in% c("casecontrol","cluster","subgroups"))
  
  # rank pathways by chosen comparison (desc), tie-break by the other two
  ord_tab <- long_df |>
    dplyr::select(pathway, contrast, neglog10padj) |>
    tidyr::pivot_wider(names_from = contrast, values_from = neglog10padj, values_fill = -Inf)
  
  pri <- sort_by
  sec <- setdiff(c("casecontrol","cluster","subgroups"), pri)
  ord_tab <- dplyr::arrange(ord_tab, dplyr::desc(.data[[pri]]),
                            dplyr::desc(.data[[sec[1]]]),
                            dplyr::desc(.data[[sec[2]]]))
  
  keep <- head(ord_tab$pathway, n = min(top_n, nrow(ord_tab)))
  long_df <- dplyr::filter(long_df, pathway %in% keep)
  long_df$pathway <- factor(long_df$pathway, levels = rev(keep))
  long_df$fill_key <- factor(long_df$fill_key, levels = .fill_levels)
  
  # within-pathway order: casecontrol (top), cluster, subgroups
  dodge_rev <- ggplot2::position_dodge2(width = 0.7, preserve = "single", reverse = TRUE)
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(y = pathway, x = neglog10padj, fill = fill_key)) +
    ggplot2::geom_col(position = dodge_rev, width = 0.65) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = .fill_vals, labels = .fill_labs, breaks = .fill_levels, name = "Direction") +
    ggplot2::labs(y = NULL, x = expression(-log[10](padj)), title = title_lab,
                  subtitle = paste0("sorted by ", sort_by, " (decreasing)")) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",    # << keep for guides='collect'
      axis.text.y = ggplot2::element_text(size = 8, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 11, margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey30", margin = ggplot2::margin(b = 6))
    )
  
  if (!is.null(xlim_max) && is.finite(xlim_max) && xlim_max > 0) {
    p <- p + ggplot2::coord_cartesian(xlim = c(0, xlim_max), expand = FALSE)
  }
  p
}

# ---- stack multiple clusters into a single figure & save
plot_multi_clusters_bars <- function(wide_nes, assignments,
                                     clusters,
                                     sort_by   = "casecontrol",
                                     top_n     = 20,
                                     out_file  = file.path("Results_subgroup_DEGs","pathways","top_pathways_per_cluster","multi_clusters.pdf"),
                                     width     = 8,
                                     height    = NULL,
                                     gap_rel   = 0.12) {
  
  # normalize sort_by length
  if (length(sort_by) == 1L) sort_by <- rep(sort_by, length(clusters))
  stopifnot(length(sort_by) == length(clusters))
  
  # shared x-axis range across selected clusters
  padj_cols <- grep("^padj_", names(wide_nes), value = TRUE)
  stopifnot(length(padj_cols) > 0)
  eps <- .Machine$double.xmin
  xlim_max <- max(-log10(pmax(unlist(wide_nes[padj_cols]), eps)), na.rm = TRUE)
  xlim_max <- ceiling(xlim_max * 1.05 * 10) / 10
  
  # build the panels
  panels <- Map(function(cl, sb) {
    .make_cluster_panel(wide_nes, assignments, cluster = cl, sort_by = sb,
                        top_n = top_n, xlim_max = xlim_max)
  }, clusters, sort_by)
  
  # interleave spacers between panels so we can control vertical gap reliably
  interleaved <- list()
  heights_vec <- c()
  for (i in seq_along(panels)) {
    interleaved[[length(interleaved) + 1]] <- panels[[i]]
    heights_vec <- c(heights_vec, 1)
    if (i < length(panels)) {
      interleaved[[length(interleaved) + 1]] <- patchwork::plot_spacer()
      heights_vec <- c(heights_vec, gap_rel) 
    }
  }
  
  combined <- patchwork::wrap_plots(interleaved, ncol = 1, heights = heights_vec) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom",
                   plot.margin = ggplot2::margin(6, 8, 6, 8))
  
  if (is.null(height)) {
    height <- length(panels) * 2.6 + (length(panels) - 1) * (2.6 * gap_rel) + 0.7
  }
  
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_file, combined, width = width, height = height, device = cairo_pdf)
  message("[save] ", normalizePath(out_file))
  
  invisible(combined)
}

# Create barplots
plot_multi_clusters_bars(
  wide_nes, assignments,
  clusters = c("vis_clus_3", "vis_clus_12", "vis_clus_6", "vis_clus_7"),
  sort_by  = c("cluster", "cluster", "casecontrol", "casecontrol"),
  top_n    = 5,
  out_file = file.path("Results_subgroup_DEGs","pathways","top_pathways_per_cluster","multi_mixed_sorting.pdf"),
  width    = 8
)




## ==== 6) Leading-edge csv for top pathways per contrast ====
lead_export <- all_fg %>%
  group_by(contrast) %>%
  slice_min(padj, n = 20, with_ties = FALSE) %>%
  transmute(contrast, pathway, NES, padj,
            leading_edge_genes = vapply(leadingEdge, function(v) paste(v, collapse = ";"), "")) %>%
  ungroup()

readr::write_csv(lead_export, file.path(PATH_DIR, "Mic_P2RY12_leading_edges_top_pathways.csv"))
message("[export] Leading edges: ", normalizePath(file.path(PATH_DIR, "Mic_P2RY12_leading_edges_top_pathways.csv")))

