library(ggplot2)
library(ggrepel)
library(gridExtra)
library(purrr)
library(reshape2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(readr)


#' Reorder a matrix based on TSP solver
#'
#' @param mat a matrix
#' @param rows Whether to cluster rows or cols
#' @param method Method passed to base::dist
#' @return A matrix sorted by rows
#' @export
order.tsp <- function(mat, rows=TRUE, method="euclidean") {
  if (rows) {
    tsp = seriation::seriate(dist(mat, method=method))
  } else {
    tsp = seriation::seriate(dist(t(mat), method=method))
  }
  ord = seriation::get_order(tsp, 1)
  if (rows) {
    return(mat[ord,])
  } else {
    return(mat[,ord])
  }
}

#' Diagonally cluster the matrix
#'
#' This function clusters a matrix
#'
#' @param mat Input matrix
#' @param ratio Ratio
#' @param cutoff Cutoff for value being counted
#' @return A list of the sorted matrix and associated columns
#' @export
diag.mat3 <- function(mat, ratio=0.5, cutoff=0.25, rows=TRUE) {
  prev.cutoff = attr(mat, "cutoff")
  if (rows) {
    ord = order(rowSums(mat > cutoff)  > ratio * ncol(mat),
                apply(mat,1,which.max), decreasing=TRUE)
    mat = mat[ord,]
    cto = apply(mat, 1, which.max)
    idx = rowSums(mat > cutoff)  > ratio * ncol(mat)
    cto[idx] = 0
  } else {
    ord = order(colSums(mat > cutoff)  > ratio * nrow(mat),
                apply(mat,2,which.max), decreasing=TRUE)
    mat = mat[,ord]
    cto = apply(mat, 2, which.max)
    idx = colSums(mat > cutoff)  > ratio * nrow(mat)
    cto[idx] = 0
    
  }
  attr(mat, "cutoff") = prev.cutoff
  if (is.null(attr(mat, "cutoff"))) {
    if (rows) {
      attr(mat, "cutoff") = list(row=cto, col=NULL)
    } else {
      attr(mat, "cutoff") = list(row=NULL, col=cto)
    }
  } else {
    if (rows) {
      attr(mat, "cutoff")$row = cto
    } else {
      attr(mat, "cutoff")$col = cto
    }
  }
  return(mat)
}

htSortMatrix <- function(M, method="euclidean", ratio=0.5, cutoff=0.25, sort=c(1,2)) {
  if (is.logical(sort)) {
    if (sort) {
      sort = c(1,2)
    } else {
      sort = c()
    }
  }
  rows = 1 %in% sort
  columns = 2 %in% sort
  if (rows) {
    M = order.tsp(M, rows=TRUE, method=method)
  }
  if (columns) {
    M = order.tsp(M, rows=FALSE, method=method)
  }
  if (rows) {
    M = diag.mat3(M, rows=TRUE, ratio=ratio, cutoff=cutoff)
  }
  if (columns) {
    M = diag.mat3(M, rows=FALSE, ratio=ratio, cutoff=cutoff)
  }
  return(M)
}



# -- tidy 
deseq_to_df_volcano <- function(res) {
  df <- as.data.frame(res)
  gene <- if ("geneName" %in% names(df)) df$geneName else rownames(df)
  # choose FDR source
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
                               fdr_thr      = 0.01,
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
  d$sig_dir[is_sig & d$log2FC > 0 & abs(d$log2FC) >= lfc_thr] <- "pos_sig"
  d$sig_dir[is_sig & d$log2FC < 0 & abs(d$log2FC) >= lfc_thr] <- "neg_sig"
  pal <- c(pos_sig = colors$pos, neg_sig = colors$neg, ns = colors$ns)
  
  p <- ggplot2::ggplot(d, aes(x = log2FC, y = mlog10FDR)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", linewidth = 0.5, color = "grey40") +
    ggplot2::geom_vline(xintercept = lfc_thr,               linetype = "dashed",  linewidth = 0.5, color = "grey40") +
    ggplot2::geom_vline(xintercept = -lfc_thr,               linetype = "dashed",  linewidth = 0.5, color = "grey40") +
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
  message(sprintf("[volcano] %s — FDR<%.2f, |logFC|>=%.1f: up=%d, down=%d; labeled up=%d, down=%d",
                  title, fdr_thr, lfc_thr, n_up, n_down, nrow(up_df), nrow(dn_df)))
  p
}

get_counts <- function(res, pAdjThreshold = 0.05) {
  df <- data.frame(
    up = sum(res$FDR < pAdjThreshold & res$log2FC > 1, na.rm = TRUE),
    down = sum(res$FDR < pAdjThreshold & res$log2FC < 1, na.rm = TRUE)
  )
  return(df)
}

plot_effect_size_comparison <- function(df1, df2,
                                        label = TRUE,
                                        n_labels = 10,
                                        colors = list(sig = "#ec2426", ns = "grey70"),
                                        fdr_thr1 = 0.01,
                                        fdr_thr2 = 0.01,
                                        lfc_thr = 1.0,
                                        point_size = 1.6,
                                        point_alpha = 0.85,
                                        seed = 123,
                                        title = "",
                                        xlab = "log2(FC) [Method 1]",
                                        ylab = "log2(FC) [Method 2]") {
  # Merge on gene
  merged <- merge(df1, df2, by = "gene", suffixes = c("_1", "_2"))
  
  # Significance
  merged$sig <- "ns"
  merged$sig[merged$FDR_1 < fdr_thr1 & merged$FDR_2 < fdr_thr2] <- "sig_both"
  
  # Scatter plot
  p <- ggplot2::ggplot(merged, aes(x = log2FC_1, y = log2FC_2)) +
    ggplot2::geom_point(aes(color = sig), size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_hline(yintercept = 0, color = "grey40", linetype = "solid") +
    ggplot2::geom_vline(xintercept = 0, color = "grey40", linetype = "solid") +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0))
  
  if (label) {
    # pick top n by FDR
    top_genes <- merged[order(merged$FDR_1 + merged$FDR_2), ]
    top_genes <- head(top_genes, n_labels)
    top_genes$clean_gene <- gsub("\\.", " ", top_genes$gene)
    
    set.seed(seed)
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = clean_gene),
      size = 3, max.overlaps = 100, box.padding = 0.35, point.padding = 0.25,
      force = 0.9, min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6
    )
  }
  
  p
}

################################################################################
############################ ANALYSIS STARTS HERE ##############################
################################################################################

MicP2RY12_DEGs<-readRDS('ResultsNoApoeCovariate/Mic P2RY12.rds')

# ---- paths
OUT_BASE  <- "Results_subgroup_DEGs"
PLOT_DIR  <- file.path(OUT_BASE, "plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


dfCaseControl <- deseq_to_df_volcano(MicP2RY12_DEGs$casecontrol)
caseControlMic <- make_deseq_volcano(dfCaseControl,
                                     title = "",
                                     colors    = list(pos = "#ec2426", neg = "#3a54a4", ns = "grey80"),
                                     label_per_side = 10
                                     )

dfSubgroup <- deseq_to_df_volcano(MicP2RY12_DEGs$subgroups)
subgroupMic <- make_deseq_volcano(dfSubgroup,
                                 title = "",
                                 colors    = list(pos = "#ad2524", neg = "#264653", ns = "grey80"),
                                 label_per_side = 10
                                 )

pdf("pdf/microglia.pdf", height=3, width=8)
grid.arrange(caseControlMic, subgroupMic, ncol=2)
dev.off()

get_counts(dfCaseControl, pAdjThreshold = 0.01)
get_counts(dfSubgroup, pAdjThreshold = 0.01)


ad_nonAD_pathways <- read.csv("~/test/AD_vs_nonAD_fgsea.csv")
cluster_pathways <- read.csv("~/test/K3_vs_K1_fgsea.csv")

ad_nonAD_pathways <- select(ad_nonAD_pathways, c("padj", "pathway"))
cluster_pathways <- select(cluster_pathways, c("padj", "pathway"))

# merge to get pathway, pAdjAdNonAD, pAdjCluster

merged_pathways <- ad_nonAD_pathways %>%
  rename(pAdjAdNonAD = padj) %>%
  full_join(cluster_pathways %>% rename(pAdjCluster = padj), by = "pathway")

merged_pathways_sig <- merged_pathways %>%
  filter(pAdjAdNonAD < 0.05 | pAdjCluster < 0.05)

merged_pathways_sig <- merged_pathways %>%
  mutate(
    significance = case_when(
      pAdjAdNonAD < 0.05 & pAdjCluster < 0.05 ~ "Significant in both",
      pAdjAdNonAD < 0.05 & (is.na(pAdjCluster) | pAdjCluster >= 0.05) ~ "Significant in case/control",
      pAdjCluster < 0.05 & (is.na(pAdjAdNonAD) | pAdjAdNonAD >= 0.05) ~ "Significant in cluster",
      TRUE ~ "Not significant"
    )
  )

subset(merged_pathways_sig, significance == "Significant in both")

################################################################################
############################## Supp Volcano ####################################
################################################################################

make_celltype_plots <- function(DE_rds,
                                fdr_thr = 0.01,
                                lfc_thr = 1.0,
                                point_size = 0.5,
                                label_per_side = 10) {
  # Load DE results
  DE_list <- readRDS(DE_rds)
  celltype_name <- gsub("Results/|\\.rds", "", DE_rds)
  
  # Convert DESeq results to volcano-ready data frames
  dfCaseControl <- deseq_to_df_volcano(DE_list$casecontrol)
  dfSubgroup    <- deseq_to_df_volcano(DE_list$subgroups)
  dfShuf        <- deseq_to_df_volcano(DE_list$shuf)
  
  # Volcano plots
  caseControlPlot <- make_deseq_volcano(dfCaseControl,
                                        title = "",
                                        point_size = point_size,
                                        colors = list(pos = "#ec2426", neg = "#3a54a4", ns = "grey80"),
                                        label_per_side = label_per_side)
  
  subgroupPlot <- make_deseq_volcano(dfSubgroup,
                                     title = "",
                                     point_size = point_size,
                                     colors = list(pos = "#ad2524", neg = "#264653", ns = "grey80"),
                                     label_per_side = label_per_side)
  
  shufPlot <- make_deseq_volcano(dfShuf,
                                 title = "",
                                 point_size = point_size,
                                 colors = list(pos = "#f77f00", neg = "#6a994e", ns = "grey80"),
                                 label_per_side = label_per_side)
  
  # Effect size comparison
  efSizePlot <- plot_effect_size_comparison(dfCaseControl, dfSubgroup,
                                            title = "",
                                            point_size = point_size,
                                            xlab = "log2(FC) [Case / Control]",
                                            ylab = "log2(FC) [Subgroup]")
  
  list(celltype = celltype_name,
       caseControl = caseControlPlot,
       subgroup    = subgroupPlot,
       shuf        = shufPlot,
       effectSize  = efSizePlot)
}

celltype_files <- c(
  'ResultsNoApoeCovariate/Exc RELN CHD7.rds',
  'ResultsNoApoeCovariate/Exc L6 THEMIS NFIA.rds',
  'ResultsNoApoeCovariate/Inh VIP ABI3BP.rds',
  'ResultsNoApoeCovariate/OPC.rds',
  'ResultsNoApoeCovariate/Oli.rds',
  'ResultsNoApoeCovariate/Mic P2RY12.rds'
)

# Generate plots for all cell types
all_celltype_plots <- lapply(celltype_files, make_celltype_plots)

# Arrange as a grid: 6 rows × 4 columns
grid_list <- do.call(c, lapply(all_celltype_plots, function(x) {
  list(x$caseControl, x$subgroup, x$shuf, x$effectSize)
}))

grid.arrange(grobs = grid_list, ncol = 4)

png("pdf/suppVolcanoes.raster.png",
    width = 12 * 300, height = 18 * 300, res = 300)
grid.arrange(grobs = grid_list, ncol = 4)
dev.off()

################################################################################
################################ DEG Counts ####################################
################################################################################
cellOrder <- c(
  "Exc L3-5 RORB PLCH1",
  "Exc L2-3 CBLN2 LINC02306",
  "Exc L6 THEMIS NFIA",
  "Exc L5-6 IT Car3",
  "Exc L4-5 RORB GABRG1",
  "Exc L3-4 RORB CUX2",
  "Exc L5-6 RORB LINC02196",
  "Exc L4-5 RORB IL1RAPL2",
  "Exc L6 CT",
  "Exc L6b",
  "Exc RELN CHD7",
  "Exc L5-6 NP",
  "Exc NRGN",
  "Inh L1-6 LAMP5 CA13",
  "Inh RYR3 TSHZ2",
  "Inh ALCAM TRPM3",
  "Inh PVALB SULF1",
  "Inh PVALB CA8 (Chandelier)",
  "Inh PVALB HTR4",
  "Inh VIP ABI3BP",
  "Inh VIP CLSTN2",
  "Inh VIP TSHZ2",
  "Inh VIP THSD7B",
  "Inh ENOX2 SPHKAP",
  "Inh L3-5 SST MAFB",
  "Inh LAMP5 NRG1 (Rosehip)",
  "Inh SORCS1 TTN",  
  "Inh PTPRK FAM19A1",
  "Inh CUX2 MSR1",
  "Ast GRM3",
  "Ast DPP10",
  "Ast CHI3L1",
  "Oli",
  "OPC",
  "Mic P2RY12"
)

files <- list.files("ResultsNoApoeCovariate", pattern = "\\.rds$", full.names = TRUE) #one file per cell type
get_cell_type <- function(file_path) {
  file_name <- basename(file_path)
  str_remove(file_name, "\\.rds$")
}

pthreshold <- 0.01
logfcThreshold <- 1

deg_summary <- lapply(files, function(f) {
  obj <- readRDS(f)
  celltype <- get_cell_type(f)
  # Convert DESeqResults objects to standard data.frames
  cc <- as.data.frame(obj$casecontrol)
  sg <- as.data.frame(obj$subgroups)
  cc <- cc %>%
    filter(!is.na(padj)) %>%
    mutate(sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold, gene = rownames(.))
  sg <- sg %>%
    filter(!is.na(padj)) %>%
    mutate(sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold, gene = rownames(.))
  merged <- full_join(
    cc %>% select(gene, log2FoldChange_cc = log2FoldChange, sig_cc = sig),
    sg %>% select(gene, log2FoldChange_sg = log2FoldChange, sig_sg = sig),
    by = "gene"
  )
  merged <- merged %>%
    mutate(category = case_when(
      sig_cc & sig_sg & sign(log2FoldChange_cc) == sign(log2FoldChange_sg) ~ "Shared",
      sig_cc & !sig_sg & sign(log2FoldChange_cc) != sign(log2FoldChange_sg) ~ "CaseControlOnly",
      !sig_cc & sig_sg ~ "SubgroupOnly",
      TRUE ~ "NotSig"
    )) %>%
    mutate(direction = ifelse(log2FoldChange_sg > 0, "Up", "Down")) %>%
    mutate(celltype = celltype) %>%
    group_by(category, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(celltype = celltype)
  merged
}) %>% bind_rows()

plot_data <- deg_summary %>%
  filter(category != "NotSig", !is.na(direction)) %>%
  mutate(
    signed_count = ifelse(direction == "Up", count, -count),
    combo = paste(category, direction, sep = "_"),
    celltype = factor(celltype, levels = rev(cellOrder))
  )

combo_levels <- c(
  "SubgroupOnly_Down", 
  "Shared_Down", 
  "CaseControlOnly_Down", 
  "SubgroupOnly_Up",
  "Shared_Up",
  "CaseControlOnly_Up"
)

plot_data <- plot_data %>%
  mutate(
    combo = factor(combo, levels = combo_levels)
  )

category_colors <- c(
  "CaseControlOnly_Up" = "#e9c46a",
  "CaseControlOnly_Down" = "#8ab17d",
  "SubgroupOnly_Up" = "#f4a261",
  "SubgroupOnly_Down" = "#2a9d8f",
  "Shared_Up" = "#e76f51",
  "Shared_Down" = "#264653"
)

make_totals <- function(df) {
  df %>%
    mutate(direction = if_else(signed_count > 0, "Up", "Down")) %>%
    group_by(celltype, direction) %>%
    summarise(total = sum(signed_count), .groups = "drop")
}
make_totals(plot_data)

sg_plot_data <- plot_data %>%
  filter(combo %in% c("Shared_Up", "Shared_Down", "SubgroupOnly_Up", "SubgroupOnly_Down"))

pdf("pdf/degSummary.main.pdf", height=5, width=10)
ggplot(sg_plot_data, aes(x = celltype, y = signed_count, fill = combo)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = category_colors, name = "Category × Direction") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of DEGs",
    title = "Case-control and subgroup DEGs within cell types"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )
dev.off()

################################################################################
############################### Combined DF ####################################
################################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)

files <- list.files("ResultsNoApoeCovariate", pattern = "\\.rds$", full.names = TRUE)
pthreshold <- 0.01
logfcThreshold <- 1

get_cell_type <- function(file_path) {
  basename(file_path) %>% str_remove("\\.rds$")
}

process_deseq_contrast <- function(df, contrast_name, pthreshold, logfcThreshold) {
  df <- as.data.frame(df) %>%
    rownames_to_column("gene") %>%  # Add gene column first!
    filter(!is.na(padj)) %>%
    mutate(
      sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold,
      direction = case_when(
        sig & log2FoldChange > 0 ~ "Up",
        sig & log2FoldChange < 0 ~ "Down",
        TRUE ~ "NotSig"
      ),
      contrast = contrast_name
    ) %>%
    select(gene, log2FoldChange, padj, sig, direction, contrast)
  
  df
}

deg_all <- map_dfr(files, function(f) {
  celltype <- get_cell_type(f)
  obj <- readRDS(f)
  map_dfr(names(obj), function(contrast) {
    df <- obj[[contrast]]
    process_deseq_contrast(df, contrast, pthreshold, logfcThreshold) %>%
      mutate(celltype = celltype)
  })
})

deg_all %>%
  mutate(
    gene = as.character(gene),
    contrast = as.character(contrast),
    celltype = as.character(celltype)
  )
deg_all_sig <- filter(deg_all, padj < 0.01 & abs(log2FoldChange) > 1)

# writing to disk for supp materials
write_tsv(deg_all_sig, file = "TableS4.AllSignificantDEGs.tsv.gz")
write_tsv(deg_all, file = "DataS1.CompleteDegAnalysisResults.tsv.gz")

# for querying individual genes
deg_all %>% filter(str_detect(gene, "SLC22")) %>% filter(celltype == "Mic P2RY12") %>% filter(contrast == "subgroups")

################################################################################
############################ Summary Statistics ################################
################################################################################

pthreshold <- 0.01
logfcThreshold <- 1

deg_stats <- lapply(files, function(f) {
  obj <- readRDS(f)
  celltype <- get_cell_type(f)
  
  # Convert to data frames and add gene names
  cc <- as.data.frame(obj$casecontrol) %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold,
      gene = rownames(.),
      direction = ifelse(log2FoldChange > 0, "Up", "Down"),
      celltype = celltype
    ) %>%
    filter(sig) %>%
    select(gene, direction, celltype)
  
  sg <- as.data.frame(obj$subgroups) %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold,
      gene = rownames(.),
      direction = ifelse(log2FoldChange > 0, "Up", "Down"),
      celltype = celltype
    ) %>%
    filter(sig) %>%
    select(gene, direction, celltype)
  
  list(casecontrol = cc, subgroup = sg)
}) 

# Extract case-control and subgroup DEGs
cc_degs <- bind_rows(lapply(deg_stats, function(x) x$casecontrol))
sg_degs <- bind_rows(lapply(deg_stats, function(x) x$subgroup))

# Summary statistics
summary_stats <- list(
  # Case-Control
  cc_total_events = nrow(cc_degs),
  cc_up_events = sum(cc_degs$direction == "Up"),
  cc_down_events = sum(cc_degs$direction == "Down"),
  cc_distinct_genes = length(unique(cc_degs$gene)),
  cc_distinct_up = length(unique(cc_degs$gene[cc_degs$direction == "Up"])),
  cc_distinct_down = length(unique(cc_degs$gene[cc_degs$direction == "Down"])),
  
  # Subgroup
  sg_total_events = nrow(sg_degs),
  sg_up_events = sum(sg_degs$direction == "Up"),
  sg_down_events = sum(sg_degs$direction == "Down"),
  sg_distinct_genes = length(unique(sg_degs$gene)),
  sg_distinct_up = length(unique(sg_degs$gene[sg_degs$direction == "Up"])),
  sg_distinct_down = length(unique(sg_degs$gene[sg_degs$direction == "Down"]))
)

# Print in a nice format
cat("=== DEG Summary Statistics ===\n\n")
cat("CASE-CONTROL COMPARISON:\n")
cat(sprintf("  Total DEG events: %d (%d up, %d down)\n", 
            summary_stats$cc_total_events, 
            summary_stats$cc_up_events, 
            summary_stats$cc_down_events))
cat(sprintf("  Distinct genes: %d (%d up, %d down)\n", 
            summary_stats$cc_distinct_genes,
            summary_stats$cc_distinct_up,
            summary_stats$cc_distinct_down))
cat(sprintf("  Avg events per distinct gene: %.2f\n", 
            summary_stats$cc_total_events / summary_stats$cc_distinct_genes))

cat("\nSUBGROUP COMPARISON:\n")
cat(sprintf("  Total DEG events: %d (%d up, %d down)\n", 
            summary_stats$sg_total_events, 
            summary_stats$sg_up_events, 
            summary_stats$sg_down_events))
cat(sprintf("  Distinct genes: %d (%d up, %d down)\n", 
            summary_stats$sg_distinct_genes,
            summary_stats$sg_distinct_up,
            summary_stats$sg_distinct_down))
cat(sprintf("  Avg events per distinct gene: %.2f\n", 
            summary_stats$sg_total_events / summary_stats$sg_distinct_genes))

# Optional: Create a data frame for easy viewing/export
summary_df <- data.frame(
  Comparison = c("Case-Control", "Subgroup"),
  Total_Events = c(summary_stats$cc_total_events, summary_stats$sg_total_events),
  Up_Events = c(summary_stats$cc_up_events, summary_stats$sg_up_events),
  Down_Events = c(summary_stats$cc_down_events, summary_stats$sg_down_events),
  Distinct_Genes = c(summary_stats$cc_distinct_genes, summary_stats$sg_distinct_genes),
  Distinct_Up = c(summary_stats$cc_distinct_up, summary_stats$sg_distinct_up),
  Distinct_Down = c(summary_stats$cc_distinct_down, summary_stats$sg_distinct_down)
)

print(summary_df)

################################################################################
########################## CPP Tx and C/C Separate #############################
################################################################################

cc_plot_data <- deg_summary %>%
  filter(category %in% c("CaseControlOnly", "Shared")) %>%
  group_by(celltype, direction) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    signed_count = ifelse(direction == "Up", count, -count),
    celltype = factor(celltype, levels = rev(cellOrder))
  )

sg_plot_data_simple <- deg_summary %>%
  filter(category %in% c("SubgroupOnly", "Shared")) %>%
  group_by(celltype, direction) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    signed_count = ifelse(direction == "Up", count, -count),
    celltype = factor(celltype, levels = rev(cellOrder))
  )


direction_colors <- c(
  "Up" = "#e76f51",
  "Down" = "#264653"
)

# Case/Control plot
p_cc <- ggplot(cc_plot_data, aes(x = celltype, y = signed_count, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = direction_colors, name = "Direction") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of DEGs",
    title = "Case/Control"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )


pdf("pdf/degSummary.main.pdf", height=5, width=7)
ggplot(sg_plot_data_simple, aes(x = celltype, y = signed_count, fill = direction)) +
  geom_col() +
  geom_text(
    data = sg_plot_data_simple %>% filter(direction == "Down"),
    aes(label = count, y = signed_count),
    hjust = 1.1,
    size = 3,
    color = "black"
  ) +
  geom_text(
    data = sg_plot_data_simple %>% filter(direction == "Up"),
    aes(label = count, y = signed_count),
    hjust = -0.1,
    size = 3,
    color = "black"
  ) +
  coord_flip() +
  scale_fill_manual(values = direction_colors, name = "Direction") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of DEGs",
    title = ""
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "none"
  )
dev.off()

################################################################################
############################## APOE Comparison #################################
################################################################################
filesNoApoeCovariate <- list.files("ResultsNoApoeCovariate", pattern = "\\.rds$", full.names = TRUE)
filesApoeCovariate <- list.files("Results", pattern = "\\.rds$", full.names = TRUE)

get_cell_type <- function(file_path) {
  file_name <- basename(file_path)
  str_remove(file_name, "\\.rds$")
}

pthreshold <- 0.01
logfcThreshold <- 1


filesNoApoeCovariate <- sort(filesNoApoeCovariate)
filesApoeCovariate <- sort(filesApoeCovariate)

# Verify matching file names
stopifnot(length(filesNoApoeCovariate) == length(filesApoeCovariate))
stopifnot(all(basename(filesNoApoeCovariate) == basename(filesApoeCovariate)))

deg_summary <- lapply(seq_along(filesNoApoeCovariate), function(i) {
  f_noApoe <- filesNoApoeCovariate[i]
  f_Apoe <- filesApoeCovariate[i]
  
  celltype <- get_cell_type(f_noApoe)
  
  # Read both files
  obj_noApoe <- readRDS(f_noApoe)
  obj_Apoe <- readRDS(f_Apoe)
  
  # Extract subgroups from both
  sg_noApoe <- as.data.frame(obj_noApoe$subgroups)
  sg_Apoe <- as.data.frame(obj_Apoe$subgroups)
  
  # Process no APOE covariate results
  sg_noApoe <- sg_noApoe %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold,
      gene = rownames(.)
    )
  
  # Process with APOE covariate results
  sg_Apoe <- sg_Apoe %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold,
      gene = rownames(.)
    )
  
  # Merge the two result sets
  merged <- full_join(
    sg_noApoe %>% select(gene, log2FoldChange_noApoe = log2FoldChange, sig_noApoe = sig),
    sg_Apoe %>% select(gene, log2FoldChange_Apoe = log2FoldChange, sig_Apoe = sig),
    by = "gene"
  )
  
  merged <- merged %>%
    mutate(category = case_when(
      sig_noApoe & sig_Apoe & sign(log2FoldChange_noApoe) == sign(log2FoldChange_Apoe) ~ "Shared",
      sig_noApoe & !sig_Apoe ~ "APOE_mediated",  # lost when controlling for APOE
      !sig_noApoe & sig_Apoe ~ "APOE_independent",  # revealed when controlling for APOE
      TRUE ~ "NotSig"
    )) %>%
    mutate(direction = ifelse(log2FoldChange_noApoe > 0, "Up", "Down")) %>%
    mutate(celltype = celltype) %>%
    group_by(category, direction, celltype) %>%
    summarise(count = n(), .groups = "drop")
  
}) %>% bind_rows()

# Prepare plot data
plot_data <- deg_summary %>%
  filter(category != "NotSig", !is.na(direction)) %>%
  mutate(
    signed_count = ifelse(direction == "Up", count, -count),
    combo = paste(category, direction, sep = "_"),
    celltype = factor(celltype, levels = rev(cellOrder))
  )

# Define combo levels
combo_levels <- c(
  "APOE_independent_Down",
  "Shared_Down",
  "APOE_mediated_Down",
  "APOE_independent_Up",
  "Shared_Up",
  "APOE_mediated_Up"
)


plot_data <- plot_data %>%
  mutate(
    combo = factor(combo, levels = combo_levels)
  )

# Define colors
category_colors <- c(
  "APOE_mediated_Up" = "#e9c46a",
  "APOE_mediated_Down" = "#8ab17d",
  "APOE_independent_Up" = "#f4a261",
  "APOE_independent_Down" = "#2a9d8f",
  "Shared_Up" = "#e76f51",
  "Shared_Down" = "#264653"
)

pdf("pdf/apoeCovariateComparison.pdf", height=5, width=10)
ggplot(plot_data, aes(x = celltype, y = signed_count, fill = combo)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = category_colors, name = "Category × Direction") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of DEGs",
    title = "Subgroup DEGs: With vs Without APOE Covariate"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )
dev.off()


################################################################################
############## Vizualizing with multiple pAdj /logFC cutoffs ###################
################################################################################

files <- list.files("ResultsNoApoeCovariate", pattern = "\\.rds$", full.names = TRUE)

get_cell_type <- function(file_path) {
  file_name <- basename(file_path)
  str_remove(file_name, "\\.rds$")
}

# Function to assign each gene to its MOST STRINGENT threshold category
assign_threshold_category <- function(f) {
  obj <- readRDS(f)
  celltype <- get_cell_type(f)
  
  sg <- as.data.frame(obj$casecontrol) %>%
    filter(!is.na(padj)) %>%
    mutate(
      gene = rownames(.),
      direction = ifelse(log2FoldChange > 0, "Up", "Down"),
      abs_lfc = abs(log2FoldChange)
    )
  
  # Assign to mutually exclusive categories (most stringent first)
  sg <- sg %>%
    mutate(
      threshold_category = case_when(
        # Most stringent: pAdj < 0.01 & |LFC| > 1
        padj < 0.01 & abs_lfc > 1 ~ "p<0.01_LFC>1",
        # pAdj < 0.01 & |LFC| > 0.5 (but not >1)
        padj < 0.01 & abs_lfc > 0.5 ~ "p<0.01_LFC>0.5",
        # pAdj < 0.01 & |LFC| > 0 (but not >0.5)
        padj < 0.01 & abs_lfc > 0 ~ "p<0.01_LFC>0",
        # pAdj < 0.05 & |LFC| > 1 (but didn't pass p<0.01)
        padj < 0.05 & abs_lfc > 1 ~ "p<0.05_LFC>1",
        # pAdj < 0.05 & |LFC| > 0.5 (but not >1)
        padj < 0.05 & abs_lfc > 0.5 ~ "p<0.05_LFC>0.5",
        # pAdj < 0.05 & |LFC| > 0 (but not >0.5)
        padj < 0.05 & abs_lfc > 0 ~ "p<0.05_LFC>0",
        # Not significant
        TRUE ~ "NotSig"
      )
    )
  
  # Count by category and direction
  counts <- sg %>%
    filter(threshold_category != "NotSig") %>%
    group_by(celltype = celltype, threshold_category, direction) %>%
    summarise(count = n(), .groups = "drop")
  
  return(counts)
}

# Run for all cell types
all_deg_counts <- lapply(files, assign_threshold_category) %>% bind_rows()

# Create combined category for plotting
plot_data <- all_deg_counts %>%
  mutate(
    signed_count = ifelse(direction == "Up", count, -count),
    combo = paste(threshold_category, direction, sep = "_"),
    celltype = factor(celltype, levels = rev(cellOrder))
  )

# Define order for categories (most stringent to least)
combo_levels <- c(
  "p<0.01_LFC>1_Down",
  "p<0.01_LFC>0.5_Down",
  "p<0.01_LFC>0_Down",
  "p<0.05_LFC>1_Down",
  "p<0.05_LFC>0.5_Down",
  "p<0.05_LFC>0_Down",
  "p<0.01_LFC>1_Up",
  "p<0.01_LFC>0.5_Up",
  "p<0.01_LFC>0_Up",
  "p<0.05_LFC>1_Up",
  "p<0.05_LFC>0.5_Up",
  "p<0.05_LFC>0_Up"
)

plot_data <- plot_data %>%
  mutate(combo = factor(combo, levels = combo_levels))

# Define 12 colors (gradient from dark to light for each direction)
# Down: blues/greens (dark to light as stringency decreases)
# Up: reds/oranges (dark to light as stringency decreases)
category_colors <- c(
  # Down-regulated (most to least stringent)
  "p<0.01_LFC>1_Down" = "#264653",      # darkest
  "p<0.01_LFC>0.5_Down" = "#2a9d8f",
  "p<0.01_LFC>0_Down" = "#8ab17d",
  "p<0.05_LFC>1_Down" = "#a8c5b8",
  "p<0.05_LFC>0.5_Down" = "#c9ddd4",
  "p<0.05_LFC>0_Down" = "#e0ece7",      # lightest
  # Up-regulated (most to least stringent)
  "p<0.01_LFC>1_Up" = "#9d0208",        # darkest
  "p<0.01_LFC>0.5_Up" = "#d00000",
  "p<0.01_LFC>0_Up" = "#e85d04",
  "p<0.05_LFC>1_Up" = "#f48c06",
  "p<0.05_LFC>0.5_Up" = "#faa307",
  "p<0.05_LFC>0_Up" = "#ffba08"         # lightest
)


ccThresholdingPlot <- ggplot(plot_data, aes(x = celltype, y = signed_count, fill = combo)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = category_colors, 
    name = "Threshold Category",
    labels = c(
      "p<0.01, |LFC|>1 (Down)",
      "p<0.01, |LFC|>0.5 (Down)",
      "p<0.01, |LFC|>0 (Down)",
      "p<0.05, |LFC|>1 (Down)",
      "p<0.05, |LFC|>0.5 (Down)",
      "p<0.05, |LFC|>0 (Down)",
      "p<0.01, |LFC|>1 (Up)",
      "p<0.01, |LFC|>0.5 (Up)",
      "p<0.01, |LFC|>0 (Up)",
      "p<0.05, |LFC|>1 (Up)",
      "p<0.05, |LFC|>0.5 (Up)",
      "p<0.05, |LFC|>0 (Up)"
    )
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of DEGs",
    title = "CPP Case/Control DEGs by Threshold Stringency"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8)
  )


### Now for subgroups
assign_threshold_category <- function(f) {
  obj <- readRDS(f)
  celltype <- get_cell_type(f)
  
  sg <- as.data.frame(obj$subgroups) %>%
    filter(!is.na(padj)) %>%
    mutate(
      gene = rownames(.),
      direction = ifelse(log2FoldChange > 0, "Up", "Down"),
      abs_lfc = abs(log2FoldChange)
    )
  
  # Assign to mutually exclusive categories (most stringent first)
  sg <- sg %>%
    mutate(
      threshold_category = case_when(
        # Most stringent: pAdj < 0.01 & |LFC| > 1
        padj < 0.01 & abs_lfc > 1 ~ "p<0.01_LFC>1",
        # pAdj < 0.01 & |LFC| > 0.5 (but not >1)
        padj < 0.01 & abs_lfc > 0.5 ~ "p<0.01_LFC>0.5",
        # pAdj < 0.01 & |LFC| > 0 (but not >0.5)
        padj < 0.01 & abs_lfc > 0 ~ "p<0.01_LFC>0",
        # pAdj < 0.05 & |LFC| > 1 (but didn't pass p<0.01)
        padj < 0.05 & abs_lfc > 1 ~ "p<0.05_LFC>1",
        # pAdj < 0.05 & |LFC| > 0.5 (but not >1)
        padj < 0.05 & abs_lfc > 0.5 ~ "p<0.05_LFC>0.5",
        # pAdj < 0.05 & |LFC| > 0 (but not >0.5)
        padj < 0.05 & abs_lfc > 0 ~ "p<0.05_LFC>0",
        # Not significant
        TRUE ~ "NotSig"
      )
    )
  
  # Count by category and direction
  counts <- sg %>%
    filter(threshold_category != "NotSig") %>%
    group_by(celltype = celltype, threshold_category, direction) %>%
    summarise(count = n(), .groups = "drop")
  
  return(counts)
}

# Run for all cell types
all_deg_counts <- lapply(files, assign_threshold_category) %>% bind_rows()

# Create combined category for plotting
plot_data <- all_deg_counts %>%
  mutate(
    signed_count = ifelse(direction == "Up", count, -count),
    combo = paste(threshold_category, direction, sep = "_"),
    celltype = factor(celltype, levels = rev(cellOrder))
  )

# Define order for categories (most stringent to least)
combo_levels <- c(
  "p<0.01_LFC>1_Down",
  "p<0.01_LFC>0.5_Down",
  "p<0.01_LFC>0_Down",
  "p<0.05_LFC>1_Down",
  "p<0.05_LFC>0.5_Down",
  "p<0.05_LFC>0_Down",
  "p<0.01_LFC>1_Up",
  "p<0.01_LFC>0.5_Up",
  "p<0.01_LFC>0_Up",
  "p<0.05_LFC>1_Up",
  "p<0.05_LFC>0.5_Up",
  "p<0.05_LFC>0_Up"
)

plot_data <- plot_data %>%
  mutate(combo = factor(combo, levels = combo_levels))

subgroupThresholdingPlot <- ggplot(plot_data, aes(x = celltype, y = signed_count, fill = combo)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = category_colors, 
    name = "Threshold Category",
    labels = c(
      "p<0.01, |LFC|>1 (Down)",
      "p<0.01, |LFC|>0.5 (Down)",
      "p<0.01, |LFC|>0 (Down)",
      "p<0.05, |LFC|>1 (Down)",
      "p<0.05, |LFC|>0.5 (Down)",
      "p<0.05, |LFC|>0 (Down)",
      "p<0.01, |LFC|>1 (Up)",
      "p<0.01, |LFC|>0.5 (Up)",
      "p<0.01, |LFC|>0 (Up)",
      "p<0.05, |LFC|>1 (Up)",
      "p<0.05, |LFC|>0.5 (Up)",
      "p<0.05, |LFC|>0 (Up)"
    )
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Number of DEGs",
    title = "CPP Subgroup DEGs by Threshold Stringency"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8)
  )

pdf("pdf/degThresholdingCounts.pdf", height=10, width=9)
grid.arrange(ccThresholdingPlot, subgroupThresholdingPlot)
dev.off()

### Summary statistics
threshold_summary <- plot_data %>%
  group_by(threshold_category) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  mutate(proportion = total_count / sum(total_count))

print(threshold_summary)


################################################################################
############################# Batch Comparison #################################
################################################################################
filesWithBatch <- list.files("ResultsNoApoeCovariate", pattern = "\\.rds$", full.names = TRUE)
filesNoBatch <- list.files("ResultsNoSeqBatchNoApoeCovariate", pattern = "\\.rds$", full.names = TRUE)

get_cell_type <- function(file_path) {
  file_name <- basename(file_path)
  str_remove(file_name, "\\.rds$")
}

filesWithBatch <- sort(filesWithBatch)
filesNoBatch <- sort(filesNoBatch)

# Verify matching file names
stopifnot(length(filesWithBatch) == length(filesNoBatch))
stopifnot(all(basename(filesWithBatch) == basename(filesNoBatch)))

compare_batch_effect <- function(contrast_name = "subgroups", pthreshold = 0.01, logfcThreshold = 1) {
  
  deg_summary <- lapply(seq_along(filesWithBatch), function(i) {
    f_withBatch <- filesWithBatch[i]
    f_noBatch <- filesNoBatch[i]
    
    celltype <- get_cell_type(f_withBatch)
    
    # Read both files
    obj_withBatch <- readRDS(f_withBatch)
    obj_noBatch <- readRDS(f_noBatch)
    
    # Extract the specified contrast from both
    res_withBatch <- as.data.frame(obj_withBatch[[contrast_name]])
    res_noBatch <- as.data.frame(obj_noBatch[[contrast_name]])
    
    # Process WITH batch covariate results
    res_withBatch <- res_withBatch %>%
      filter(!is.na(padj)) %>%
      mutate(
        sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold,
        gene = rownames(.)
      )
    
    # Process WITHOUT batch covariate results
    res_noBatch <- res_noBatch %>%
      filter(!is.na(padj)) %>%
      mutate(
        sig = padj < pthreshold & abs(log2FoldChange) > logfcThreshold,
        gene = rownames(.)
      )
    
    # Merge the two result sets
    merged <- full_join(
      res_withBatch %>% select(gene, log2FoldChange_withBatch = log2FoldChange, sig_withBatch = sig),
      res_noBatch %>% select(gene, log2FoldChange_noBatch = log2FoldChange, sig_noBatch = sig),
      by = "gene"
    )
    
    merged <- merged %>%
      mutate(category = case_when(
        sig_withBatch & sig_noBatch & sign(log2FoldChange_withBatch) == sign(log2FoldChange_noBatch) ~ "Shared",
        !sig_withBatch & sig_noBatch ~ "Batch_confounded",  # lost when controlling for batch (likely false positives)
        sig_withBatch & !sig_noBatch ~ "Batch_independent",  # revealed when controlling for batch
        TRUE ~ "NotSig"
      )) %>%
      mutate(direction = ifelse(log2FoldChange_withBatch > 0, "Up", "Down")) %>%
      mutate(celltype = celltype) %>%
      group_by(category, direction, celltype) %>%
      summarise(count = n(), .groups = "drop")
    
  }) %>% bind_rows()
  
  return(deg_summary)
}

###### Plotting function ######
plot_batch_comparison <- function(deg_summary, contrast_label, pthreshold, logfcThreshold) {
  
  # Prepare plot data
  plot_data <- deg_summary %>%
    filter(category != "NotSig", !is.na(direction)) %>%
    mutate(
      signed_count = ifelse(direction == "Up", count, -count),
      combo = paste(category, direction, sep = "_"),
      celltype = factor(celltype, levels = rev(cellOrder))
    )
  
  # Define combo levels
  combo_levels <- c(
    "Batch_independent_Down",
    "Shared_Down",
    "Batch_confounded_Down",
    "Batch_independent_Up",
    "Shared_Up",
    "Batch_confounded_Up"
  )
  
  plot_data <- plot_data %>%
    mutate(
      combo = factor(combo, levels = combo_levels)
    )
  
  # Define colors
  category_colors <- c(
    "Batch_confounded_Up" = "#e9c46a",
    "Batch_confounded_Down" = "#8ab17d",
    "Batch_independent_Up" = "#f4a261",
    "Batch_independent_Down" = "#2a9d8f",
    "Shared_Up" = "#e76f51",
    "Shared_Down" = "#264653"
  )
  
  p <- ggplot(plot_data, aes(x = celltype, y = signed_count, fill = combo)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = category_colors, name = "Category × Direction") +
    theme_classic() +
    labs(
      x = NULL,
      y = "Number of DEGs",
      title = paste0(contrast_label, " DEGs (pAdj<", pthreshold, ", |LFC|>", logfcThreshold, 
                     "): With vs Without Batch Covariate")
    ) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )
  
  return(list(plot = p, data = plot_data))
}

###### Run for both contrasts ######
pthreshold <- 0.01
logfcThreshold <- 1
deg_summary_subgroups <- compare_batch_effect("subgroups", pthreshold, logfcThreshold)
subgroups_result <- plot_batch_comparison(deg_summary_subgroups, "CPP Subgroups", pthreshold, logfcThreshold)

# summary statistics
sg_summary <- deg_summary_subgroups %>%
  filter(category != "NotSig") %>%
  group_by(category) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  mutate(proportion = total_count / sum(total_count))
print(sg_summary)

pthreshold <- 0.05
logfcThreshold <- 0
deg_summary_casecontrol <- compare_batch_effect("casecontrol", pthreshold, logfcThreshold)
casecontrol_result <- plot_batch_comparison(deg_summary_casecontrol, "Case/Control", pthreshold, logfcThreshold)

cc_summary <- deg_summary_casecontrol %>%
  filter(category != "NotSig") %>%
  group_by(category) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  mutate(proportion = total_count / sum(total_count))
print(cc_summary)

pdf("pdf/degBatchCovariateComparison.pdf", height=10, width=9)
grid.arrange(subgroups_result$plot, casecontrol_result$plot)
dev.off()

################################################################################
################################ DEG Heatmap ###################################
################################################################################
cellOrder <- c(
  "Exc L3-5 RORB PLCH1",
  "Exc L2-3 CBLN2 LINC02306",
  "Exc L6 THEMIS NFIA",
  "Exc L5-6 IT Car3",
  "Exc L4-5 RORB GABRG1",
  "Exc L3-4 RORB CUX2",
  "Exc L5-6 RORB LINC02196",
  "Exc L4-5 RORB IL1RAPL2",
  "Exc L6 CT",
  "Exc L6b",
  "Exc RELN CHD7",
  "Exc L5-6 NP",
  "Exc NRGN",
  "Inh L1-6 LAMP5 CA13",
  "Inh RYR3 TSHZ2",
  "Inh ALCAM TRPM3",
  "Inh PVALB SULF1",
  "Inh PVALB CA8 (Chandelier)",
  "Inh PVALB HTR4",
  "Inh VIP ABI3BP",
  "Inh VIP CLSTN2",
  "Inh VIP TSHZ2",
  "Inh VIP THSD7B",
  "Inh ENOX2 SPHKAP",
  "Inh L3-5 SST MAFB",
  "Inh LAMP5 NRG1 (Rosehip)",
  "Inh PTPRK FAM19A1",
  "Inh CUX2 MSR1",
  "Inh ALCAM TRPM3",
  "Ast GRM3",
  "Ast DPP10",
  "Ast CHI3L1",
  "Oli",
  "OPC",
  "Mic P2RY12"
)


files <- list.files("Results", pattern = "\\.rds$", full.names = TRUE)
get_cell_type <- function(file_path) {
  file_name <- basename(file_path)
  str_remove(file_name, "\\.rds$")
}

deg_list <- map(files, function(f) {
  deg_df <- deseq_to_df_volcano(readRDS(f)$subgroups)
  deg_df$cellType <- get_cell_type(f)
  deg_df
})

# name the list
names(deg_list) <- map_chr(files, get_cell_type)
deg_all <- bind_rows(deg_list)

# pivot to wide format: genes x cellTypes (log2FC as values)
log2FC_matrix <- deg_all %>%
  select(gene, cellType, log2FC) %>%
  pivot_wider(names_from = cellType, values_from = log2FC)

# set gene names as rownames and remove gene column
log2FC_matrix <- as.data.frame(log2FC_matrix)
rownames(log2FC_matrix) <- log2FC_matrix$gene
log2FC_matrix$gene <- NULL

# optional: replace NAs with 0 (or leave as NA if you prefer)
log2FC_matrix[is.na(log2FC_matrix)] <- 0

log2FC_matrix_ordered <- log2FC_matrix[, cellOrder]

# filter to genes significant (padj < 0.05) in at least one cell type
sig_genes <- deg_all %>%
  filter(!is.na(padj) & padj < 1e-18) %>%
  pull(gene) %>%
  unique()

log2FC_sig_matrix <- log2FC_matrix[rownames(log2FC_matrix) %in% sig_genes, ]





my_colors <- colorRampPalette(c("#0077b6", "white", "#c1121f"))(100)
pheatmap(log2FC_matrix,
         color = my_colors,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         breaks = seq(-2, 2, length.out = 101)
         )



################################################################################
####################### GWAS INTEGRATION #######################################
################################################################################
# from Uffelmann et al MedRxiv
gwas_genes <- readLines("gwasList.txt")
alias_map <- c(
  "LAMPTM5" = "LAPTM5", # evidently just a misprint in the GWAS
  "GBA1" = "GBA", # alias gene name, GBA in expression, GBA1 in GWAS, refers to same gene
  "FAM105B" = "OTULIN", # again, alternative gene name
  "KIAA0125" = "FAM30A", #alt gene name
  "C6orf57" = "SDHAF4" #alt gene name, still absent gene dataset, LOST
)

gwas_genes_fixed <- ifelse(gwas_genes %in% names(alias_map),
                           alias_map[gwas_genes],
                           gwas_genes)

subset_matrix <- log2FC_matrix[rownames(log2FC_matrix) %in% gwas_genes_fixed, ]

# let's get a matrix of FDR to draw significance stars:
deg_sub <- deg_all %>% filter(gene %in% rownames(subset_matrix))

# Pivot to same shape as log2FC_matrix
FDR_matrix <- deg_sub %>%
  select(gene, cellType, FDR) %>%
  pivot_wider(names_from = cellType, values_from = FDR) %>%
  column_to_rownames("gene")

sig_matrix <- FDR_matrix
sig_matrix[,] <- ""
sig_matrix[FDR_matrix <= 0.01] <- "*"
sig_matrix <- sig_matrix[rownames(subset_matrix), colnames(subset_matrix)]

my_colors <- colorRampPalette(c("#0077b6", "white", "#c1121f"))(100)
subset_matrix_ordered <- subset_matrix[, cellOrder]
pdf("pdf/gwasPanel.pdf", height=6, width=24)
pheatmap(t(subset_matrix_ordered),
         color = my_colors,
         breaks = seq(-2, 2, length.out = 101),
         display_numbers = t(sig_matrix),
         number_color = "black",
         cluster_rows = FALSE,  # cell types
         cluster_cols = TRUE)  # genes
dev.off()

################################################################################
########################## Trying with log2FC * minusLogFDR ####################
################################################################################
cellOrder <- c(
  "Exc L3-5 RORB PLCH1",
  "Exc L2-3 CBLN2 LINC02306",
  "Exc L6 THEMIS NFIA",
  "Exc L5-6 IT Car3",
  "Exc L4-5 RORB GABRG1",
  "Exc L3-4 RORB CUX2",
  "Exc L5-6 RORB LINC02196",
  "Exc L4-5 RORB IL1RAPL2",
  "Exc L6 CT",
  "Exc L6b",
  "Exc RELN CHD7",
  "Exc L5-6 NP",
  "Exc NRGN",
  "Inh L1-6 LAMP5 CA13",
  "Inh RYR3 TSHZ2",
  "Inh ALCAM TRPM3",
  "Inh PVALB SULF1",
  "Inh PVALB CA8 (Chandelier)",
  "Inh PVALB HTR4",
  "Inh VIP ABI3BP",
  "Inh VIP CLSTN2",
  "Inh VIP TSHZ2",
  "Inh VIP THSD7B",
  "Inh ENOX2 SPHKAP",
  "Inh L3-5 SST MAFB",
  "Inh LAMP5 NRG1 (Rosehip)",
  "Inh PTPRK FAM19A1",
  "Inh CUX2 MSR1",
  "Inh ALCAM TRPM3",
  "Ast GRM3",
  "Ast DPP10",
  "Ast CHI3L1",
  "Oli",
  "OPC",
  "Mic P2RY12"
)

fdr_thr <- 0.01
lfc_thr <- 1


files <- list.files("Results", pattern = "\\.rds$", full.names = TRUE)
get_cell_type <- function(file_path) {
  file_name <- basename(file_path)
  str_remove(file_name, "\\.rds$")
}

deg_list <- map(files, function(f) {
  deg_df <- deseq_to_df_volcano(readRDS(f)$subgroups)
  deg_df$cellType <- get_cell_type(f)
  deg_df
})

# name the list
names(deg_list) <- map_chr(files, get_cell_type)
deg_all <- bind_rows(deg_list)

# pivot to wide format: genes x cellTypes (log2FC as values)
log2FC_matrix <- deg_all %>%
  select(gene, cellType, log2FC) %>%
  pivot_wider(names_from = cellType, values_from = log2FC)

# set gene names as rownames and remove gene column
log2FC_matrix <- as.data.frame(log2FC_matrix)
rownames(log2FC_matrix) <- log2FC_matrix$gene
log2FC_matrix$gene <- NULL

# optional: replace NAs with 0 (or leave as NA if you prefer)
log2FC_matrix[is.na(log2FC_matrix)] <- 0

log2FC_matrix_ordered <- log2FC_matrix[, cellOrder]

gwas_df <- read.delim("gwasList.withSig.txt", header = TRUE, sep = "\t")
colnames(gwas_df) <- c("gene", "p_value", "minusLogP")

# Define gene alias mapping
alias_map <- c(
  "LAMPTM5" = "LAPTM5",
  "GBA1" = "GBA",
  "FAM105B" = "OTULIN",
  "KIAA0125" = "FAM30A",
  "C6orf57" = "SDHAF4"
)

# Fix aliases
gwas_df$gene_fixed <- ifelse(gwas_df$gene %in% names(alias_map),
                             alias_map[gwas_df$gene],
                             gwas_df$gene)

# Subset DEG data to GWAS genes
deg_filtered <- deg_all %>% filter(gene %in% gwas_df$gene_fixed)

# Pivot log2FC and compute fcXfdr
fcXfdr_matrix <- deg_filtered %>%
  filter(!is.na(log2FC) & !is.na(mlog10FDR)) %>%
  mutate(fcXfdr = log2FC * mlog10FDR) %>%
  select(gene, cellType, fcXfdr) %>%
  pivot_wider(names_from = cellType, values_from = fcXfdr) %>%
  column_to_rownames("gene")

# Subset to GWAS genes in correct order
gwas_genes_fixed <- gwas_df$gene_fixed
fcXfdr_subset <- fcXfdr_matrix[rownames(fcXfdr_matrix) %in% gwas_genes_fixed, ]
fcXfdr_subset_ordered <- fcXfdr_subset[, cellOrder]

# Prepare significance stars
FDR_matrix <- deg_filtered %>%
  select(gene, cellType, FDR) %>%
  pivot_wider(names_from = cellType, values_from = FDR) %>%
  column_to_rownames("gene")
log2FC_matrix_filtered <- deg_filtered %>%
  select(gene, cellType, log2FC) %>%
  pivot_wider(names_from = cellType, values_from = log2FC) %>%
  column_to_rownames("gene")

sig_matrix <- FDR_matrix
sig_matrix[,] <- ""
sig_matrix[FDR_matrix <= fdr_thr & abs(log2FC_matrix_filtered) >= lfc_thr] <- "*"
sig_matrix <- sig_matrix[rownames(fcXfdr_subset), colnames(fcXfdr_subset)]

# Cap GWAS minusLogP for visualization
logPCap <- 30
annotation_row <- gwas_df %>%
  filter(gene_fixed %in% rownames(fcXfdr_subset)) %>%
  mutate(minusLogP = pmin(minusLogP, logPCap)) %>%
  select(gene_fixed, minusLogP) %>%
  column_to_rownames("gene_fixed")

# Color palette
my_colors <- colorRampPalette(c("#0077b6", "white", "#c1121f"))(100)

# Plot heatmap
pheatmap(t(fcXfdr_subset_ordered),
         color = my_colors,
         breaks = seq(-4, 4, length.out = 101),
         display_numbers = t(sig_matrix),
         number_color = "black",
         cluster_rows = FALSE,  # cell types
         cluster_cols = TRUE,   # genes
         annotation_col = annotation_row)


################################################################################
########################## Subheatmaps ####################
################################################################################
sig_logical <- sig_matrix == "*"
n_sig_per_gene <- rowSums(sig_logical)
direction_score <- rowSums(fcXfdr_subset_ordered * sig_logical)
gene_order <- order(direction_score, decreasing = TRUE)
pheatmap(t(fcXfdr_subset_ordered)[, gene_order],
         color = my_colors,
         breaks = seq(-4, 4, length.out = 101),
         display_numbers = t(sig_matrix)[, gene_order],
         number_color = "black",
         cluster_rows = FALSE,  # cell types
         cluster_cols = FALSE,   # genes
         annotation_col = annotation_row)


################################################################################
################### TAKE 2 #####################################################
################################################################################

sig_logical <- sig_matrix == "*"

# Count number of  cell types up vs down
up_count <- rowSums((fcXfdr_subset_ordered > 0))
down_count <- rowSums((fcXfdr_subset_ordered < 0))

# Net direction score
net_score <- up_count - down_count

# Assign broad group
gene_group <- rep("Neutral", length(net_score))
gene_group[net_score > 0] <- "Broad Up"
gene_group[net_score < 0] <- "Broad Down"
gene_group[up_count + down_count == 0] <- "Not Significant"

table(gene_group)

gene_order <- c(
  names(sort(up_count[gene_group == "Broad Up"], decreasing = TRUE)),
  names(sort(net_score[gene_group == "Neutral"], decreasing = TRUE)),
  names(sort(down_count[gene_group == "Broad Down"], decreasing = FALSE)),
  names(which(gene_group == "Not Significant"))
)

pdf("pdf/gwasPanel.v2.pdf", height=6, width=24)
pheatmap(
  t(fcXfdr_subset_ordered)[, gene_order],
  color = my_colors,
  breaks = seq(-4, 4, length.out = 101),
  display_numbers = t(sig_matrix)[, gene_order],
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_row,
  fontsize_number = 12  # make numbers larger
)
dev.off()
