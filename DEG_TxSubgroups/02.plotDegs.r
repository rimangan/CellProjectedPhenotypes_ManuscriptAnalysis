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


deseq_to_df_volcano <- function(res) {
  df <- as.data.frame(res)
  gene <- if ("geneName" %in% names(df)) df$geneName else rownames(df)
  
  out <- data.frame(
    gene     = as.character(gene),
    log2FC   = suppressWarnings(as.numeric(df$log2FoldChange)),
    pval     = suppressWarnings(as.numeric(df$pvalue)),
    padj     = suppressWarnings(as.numeric(df$padj)),
    FDR      = suppressWarnings(as.numeric(df$padj)),
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
    # keep only genes passing the primary significance definition
    sig <- subset(side_df, FDR < fdr_thr & abs(log2FC) >= lfc_thr)
    
    # order by strongest statistical evidence, then effect size
    sig <- sig[order(sig$FDR, -abs(sig$log2FC)), , drop = FALSE]
    
    # treat n as a ceiling
    head(sig, min(n, nrow(sig)))
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
    up = sum(res$FDR < pAdjThreshold & res$log2FC >= 1, na.rm = TRUE),
    down = sum(res$FDR < pAdjThreshold & res$log2FC <= -1, na.rm = TRUE)
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
  sg <- as.data.frame(obj$subgroups) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(sig = padj < pthreshold & abs(log2FoldChange) >= logfcThreshold)
  cc <- as.data.frame(obj$casecontrol) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(sig = padj < pthreshold & abs(log2FoldChange) >= logfcThreshold)
  merged <- full_join(
    cc %>% select(gene, log2FoldChange_cc = log2FoldChange, sig_cc = sig),
    sg %>% select(gene, log2FoldChange_sg = log2FoldChange, sig_sg = sig),
    by = "gene"
  ) %>%
    mutate(
      sig_cc = coalesce(sig_cc, FALSE),
      sig_sg = coalesce(sig_sg, FALSE)
    )
  
  merged <- merged %>%
    mutate(category = case_when(
      sig_cc & sig_sg  ~ "Shared",
      sig_cc & !sig_sg ~ "CaseControlOnly",
      !sig_cc & sig_sg ~ "SubgroupOnly",
      TRUE             ~ "NotSig"
    )) %>%
    mutate(direction = case_when(
      !is.na(log2FoldChange_sg) ~ ifelse(log2FoldChange_sg > 0, "Up", "Down"),
      !is.na(log2FoldChange_cc) ~ ifelse(log2FoldChange_cc > 0, "Up", "Down"),
      TRUE                      ~ NA_character_
    )) %>%
    mutate(celltype = celltype) %>%
    group_by(category, direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(celltype = celltype)
  merged
}) %>% bind_rows()


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
############################## Supp Volcano ####################################
################################################################################

make_celltype_plots <- function(DE_rds,
                                fdr_thr = 0.01,
                                lfc_thr = 1.0,
                                point_size = 0.5,
                                label_per_side = 5) {
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
deg_all_sig <- filter(deg_all, padj < 0.01 & abs(log2FoldChange) >= 1)

# writing to disk for supp materials
write_tsv(deg_all_sig, file = "TableS4.AllSignificantDEGs.tsv.gz")
write_tsv(deg_all, file = "DataS1.CompleteDegAnalysisResults.tsv.gz")

# for querying individual genes
deg_all %>% filter(str_detect(gene, "SLC22")) %>% filter(celltype == "Mic P2RY12") %>% filter(contrast == "subgroups")




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