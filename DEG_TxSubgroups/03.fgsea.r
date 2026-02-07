library(dplyr)
library(fgsea)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(pheatmap)
library(forcats)
library(gridExtra)
library(ggpubr)


OUT_BASE   <- "Results"
PATH_DIR   <- "Pathways"


make_ranks <- function(deseq, method="stat") {
  df <- as.data.frame(deseq)
  df <- subset(df, is.finite(df$padj))
  df$log2FoldChange <- as.numeric(as.character(df$log2FoldChange))
  df$minusLogP <- -log10(df$pvalue)
  if (method == "stat") {
    rankings <- df$stat
  } else {
    rankings <- sign(df$log2FoldChange)*(df$minusLogP)
  }
  names(rankings) <- df$geneName
  rankings <- sort(rankings, decreasing = TRUE)
  return(rankings)
}

run_fgsea_on_res <- function(res, label,
                             gsets, 
                             minSize = 5,
                             maxSize = 1000,
                             cell_type = NULL,
                             method = "stat"
) {
  ranks <- make_ranks(res, method)
  fg    <- fgsea::fgsea(pathways = gsets, 
                        stats = ranks,
                        minSize = minSize,
                        maxSize = maxSize)
  fg <- dplyr::arrange(fg, padj) %>%
    dplyr::mutate(contrast = label,
                  direction = ifelse(NES >= 0, "up", "down"))
  
  # BUILD FILENAME WITH CELL_TYPE
  filename_parts <- c(
    if (!is.null(cell_type)) cell_type else NULL,
    gsub("[^A-Za-z0-9.-]", "", label),
    "fgsea.csv"
  )
  out_csv <- file.path(PATH_DIR, paste(filename_parts, collapse = "_"))
  
  readr::write_csv(fg, out_csv)
  message("[fgsea] Saved: ", normalizePath(out_csv))
  fg
}

build_gene_sets <- function() {
  msig_h <- msigdbr(species = "Homo sapiens", collection = "H") %>% # hallmarks
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()
  
  #msig_CP <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP") %>% # canonical pathways
  #  dplyr::select(gs_name, gene_symbol) %>%
  #  distinct()
  
  msig_go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP") %>% # go biological process
    dplyr::select(gs_name, gene_symbol) %>%
    distinct()
  
  gsets_h <- split(msig_h$gene_symbol, msig_h$gs_name)
  #gsets_CP <- split(msig_CP$gene_symbol, msig_CP$gs_name)
  gsets_go <- split(msig_go$gene_symbol, msig_go$gs_name)
  
  gsets <- c(gsets_h, gsets_go)
  return(gsets)
}

run_fgsea_for_cell <- function(file, method="stat") {
  cell_type <- gsub("\\.rds", "", basename(file))
  obj <- readRDS(file)
  # Run FGSEA for both contrasts
  fg_case <- run_fgsea_on_res(as.data.frame(obj$casecontrol), "CaseControl", gsets, cell_type = cell_type, method = method)
  fg_sub  <- run_fgsea_on_res(as.data.frame(obj$subgroups), "PathAD.TxAD_vs_PathNonAD.TxNonAD", gsets, cell_type = cell_type, method = method)
  # Merge, compute -log10(p), add cell type
  all_fg <- bind_rows(fg_case, fg_sub) %>%
    mutate(
      minusLogPAdj = -log10(padj),
      minusLogP = -log10(pval),
      cellType = cell_type
    )
  
  # Return both long and wide formats
  fg_wide <- all_fg %>%
    select(pathway, contrast, NES, minusLogPAdj) %>%
    pivot_wider(
      names_from = contrast,
      values_from = c(NES, minusLogPAdj),
      names_sep = "_"
    ) %>%
    mutate(cellType = cell_type)
  
  list(long = all_fg, wide = fg_wide)
}

################################################################################
############################## ANALYSIS STARTS HERE ############################
################################################################################

rds_files <- list.files("ResultsNoApoeCovariate", pattern = "\\.rds$", full.names = TRUE)
gsets <- build_gene_sets()

fg_results <- map(rds_files, run_fgsea_for_cell)
all_fg_long <- bind_rows(map(fg_results, "long"))
all_fg_wide <- bind_rows(map(fg_results, "wide"))

celltype_to_plot <- c(
  'Exc RELN CHD7',
  'Exc L6 THEMIS NFIA',
  'Inh VIP ABI3BP',
  'OPC',
  'Oli',
  'Mic P2RY12'
)

minusLogPAdjThreshold <- 4

plot_fgsea_for_cell <- function(celltype_name) {
  df <- all_fg_wide %>%
    filter(cellType == celltype_name) %>%
    subset(minusLogPAdj_CaseControl > minusLogPAdjThreshold |
             minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD > minusLogPAdjThreshold)
  
  ggplot(df, aes(
    x = minusLogPAdj_CaseControl,
    y = minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD,
    label = pathway,
    color = minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD
  )) +
    geom_point(alpha = 0.8, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey80") +
    geom_hline(yintercept = minusLogPAdjThreshold, linetype = "dashed") +
    geom_vline(xintercept = minusLogPAdjThreshold, linetype = "dashed") +
    geom_text_repel(size = 2, max.overlaps = 20,
                    aes(color = minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD)) +
    scale_color_gradient(low = "grey70", high = "maroon") +
    labs(
      x = "-log10(P) (Case vs Control)",
      y = "-log10(P) (CPP Subgroups)",
      title = celltype_name
    ) +
    theme_classic() +
    theme(legend.position = "none")
}

all_cell_types <- unique(all_fg_wide$cellType)
fgsea_plots <- map(all_cell_types, plot_fgsea_for_cell)
names(fgsea_plots) <- all_cell_types

pdf("pdf/fgsea_all_cell_types.pdf", height=40, width=30)
grid.arrange(
  grobs = fgsea_plots, ncol = 4
)
dev.off()

################################################################################
############################## Top 5 ###########################################
################################################################################

library(dplyr)
library(ggplot2)
library(forcats)
library(tidytext)  # for reorder_within()

top5_per_cell <- all_fg_wide %>%
  select(
    pathway,
    cellType,
    NES = NES_PathAD.TxAD_vs_PathNonAD.TxNonAD,
    minusLogPAdj = `minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD`
  ) %>%
  filter(is.finite(minusLogPAdj), is.finite(NES)) %>%
  group_by(cellType) %>%
  slice_max(order_by = minusLogPAdj, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pathway = reorder_within(pathway, minusLogPAdj, cellType))

ggplot(top5_per_cell, aes(x = minusLogPAdj, y = pathway, fill = NES)) +
  geom_col() +
  facet_wrap(~ cellType, scales = "free_y") +
  scale_y_reordered() +  # ensures proper within-facet ordering
  scale_fill_gradient2(
    low = "#0077b6", mid = "white", high = "#c1121f", midpoint = 0
  ) +
  geom_vline(xintercept = 2, linetype = "dashed") +
  labs(
    x = expression(-log[10]("q-value")),
    y = "Pathway",
    fill = "NES"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 6)
  )


################################################################################
##########################   Num. Pathway   ####################################
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(patchwork)
})

# ---------- Settings ----------
sig_thr  <- 2
PLOT_DIR <- if (exists("PLOT_DIR")) PLOT_DIR else "pdf"
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)




COL_CPP <- "#036666"
COL_SH  <- "#56ab91"
COL_CC  <- "#99e2b4"

COL_CPP <- "#e9c46a"
COL_SH  <- "#e76f51"
COL_CC  <- "#264653"

# Pastel enrichment colors (CC red, Shared purple, CPP blue)
COL_CC <- "#749c75"
COL_SH  <- "#b2bd7e"
COL_CPP  <- "#e9d985"


# Major-type strip colors
COL_MAJOR <- c(
  "Ex. Neuron" = "#2E7D32",
  "Inh. Neuron"= "#A5D6A7",
  "Ast"        = "#F06292",
  "OPC"        = "#FFB74D",
  "Oli"        = "#8D6E63",
  "Mic"        = "#B39DDB",
  "Other"      = "#FFE082"
)

# --- tweak knobs (spacing) ---
STRIP_X  <- 1.20   
LABEL_X  <- 0.55   
X_LIMS   <- c(-0.05, 2.10)

short_cell <- function(x){
  x <- gsub("Exc ", "E ", x, fixed = TRUE)
  x <- gsub("Inh ", "I ", x, fixed = TRUE)
  x <- gsub("Mic ", "M ",  x, fixed = TRUE)
  x
}
infer_major <- function(ct){
  if (grepl("^Exc\\b", ct))        return("Ex. Neuron")
  if (grepl("^Inh\\b", ct))        return("Inh. Neuron")
  if (grepl("^Ast\\b", ct))        return("Ast")
  if (grepl("^OPC\\b", ct))        return("OPC")
  if (grepl("^Oli\\b|^Oligo", ct)) return("Oli")
  if (grepl("^Mic\\b", ct))        return("Mic")
  "Other"
}

# ---------- Build counts once ----------
counts_input <- all_fg_wide %>%
  transmute(cellType, pathway,
            mlogq_case = minusLogPAdj_CaseControl,
            mlogq_cpp  = minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD) %>%
  mutate(sig_case = is.finite(mlogq_case) & mlogq_case >= sig_thr,
         sig_cpp  = is.finite(mlogq_cpp)  & mlogq_cpp  >= sig_thr,
         group = dplyr::case_when(
           sig_case & sig_cpp ~ "Shared",
           sig_case           ~ "C/C only",
           sig_cpp            ~ "CPP only",
           TRUE               ~ "None")) %>%
  filter(group != "None") %>%
  distinct(cellType, group, pathway) %>%
  dplyr::count(cellType, group, name = "n") %>%
  tidyr::complete(
    cellType,
    group = factor(c("C/C only","Shared","CPP only"),
                   levels = c("C/C only","Shared","CPP only")),
    fill = list(n = 0)
  ) %>%
  mutate(group = factor(group, levels = c("C/C only","Shared","CPP only")))

# ---------- Standard cell-type order (NO sorting by totals) ----------
# 1) Prefer factor levels if already set, otherwise first appearance in all_fg_wide
base_order <- if (is.factor(all_fg_wide$cellType)) {
  levels(all_fg_wide$cellType)
} else {
  unique(as.character(all_fg_wide$cellType))
}

# 2) Keep majors grouped in a standard major order, preserve within-major base order
order_levels <- tibble(cellType = base_order) %>%
  mutate(major = vapply(cellType, infer_major, character(1)),
         major_order = factor(major, levels = names(COL_MAJOR)),
         within = match(cellType, base_order)) %>%
  arrange(major_order, within) %>%
  pull(cellType)

counts_input <- counts_input %>%
  mutate(
    cellType       = factor(as.character(cellType), levels = order_levels),
    cellType_short = factor(short_cell(as.character(cellType)),
                            levels = short_cell(as.character(order_levels))),
    major          = vapply(as.character(cellType), infer_major, character(1))
  )

# Shared y-scale (reverse so first at top)
y_levels <- rev(levels(counts_input$cellType))

# ---------- Left strip (labels + tiles closer to barplots) ----------
strip_df <- counts_input %>% distinct(cellType, major) %>%
  mutate(y = factor(cellType, levels = y_levels))

lab_df <- strip_df %>%
  mutate(major_order = factor(major, levels = names(COL_MAJOR)),
         y_num = as.integer(y)) %>%
  arrange(major_order, desc(y_num)) %>%
  group_by(major) %>% summarise(ymin = min(y_num), ymax = max(y_num), .groups = "drop") %>%
  mutate(y_mid = (ymin + ymax)/2)

p_left <- ggplot(strip_df, aes(x = STRIP_X, y = factor(cellType, levels = y_levels), fill = major)) +
  geom_tile(width = 1.25, height = 0.9, color = NA) +
  scale_fill_manual(values = COL_MAJOR, guide = "none") +
  scale_y_discrete(limits = y_levels) +
  scale_x_continuous(limits = X_LIMS) +
  geom_text(data = lab_df, aes(x = LABEL_X, y = y_mid, label = major),
            hjust = 1, size = 3.0, fontface = "bold", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  theme_void() +
  # smaller right margin to reduce gap to the barplots
  theme(plot.margin = margin(10, 0, 10, 16))

# ---------- COUNTS (horizontal; WANT: CC → Shared → CPP) ----------
counts_plot <- counts_input %>% mutate(y = factor(cellType, levels = y_levels))

p_counts <- ggplot(counts_plot, aes(y = y, x = n, fill = group)) +
  geom_col(
    position = position_stack(reverse = TRUE),  # <-- left-to-right order
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c("C/C only"=COL_CC, "Shared"=COL_SH, "CPP only"=COL_CPP),
    breaks = c("C/C only","Shared","CPP only"),
    name = "Enrichment Type"
  ) +
  labs(x = "Number of Significant Pathways", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "left",
    legend.direction = "vertical",
    legend.title = element_text(size = 7.5),
    legend.text  = element_text(size = 7),
    legend.key.size = unit(3.2, "mm"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # shrink left margin so bars sit closer to the major strip
    plot.margin = margin(10, 10, 10, 0)
  ) +
  scale_y_discrete(limits = y_levels, expand = c(0,0))

# ---------- PROPORTIONS (horizontal; rows sum to 1; WANT: CC → Shared → CPP) ----------
totals <- counts_input %>% group_by(cellType) %>% summarise(total = sum(n), .groups = "drop")

props_plot <- counts_input %>%
  left_join(totals, by = "cellType") %>%
  mutate(prop = ifelse(total > 0, n / total, 0),
         y = factor(cellType, levels = y_levels))

p_props <- ggplot(props_plot, aes(y = y, x = prop, fill = group)) +
  geom_col(
    position = position_stack(reverse = TRUE),  # <-- key for left-to-right order
    linewidth = 0.25
  ) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(
    values = c("C/C only"=COL_CC, "Shared"=COL_SH, "CPP only"=COL_CPP),
    breaks = c("C/C only","Shared","CPP only"),
    guide = "none"
  ) +
  labs(x = "Proportion of significant pathways", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 10, 10, 0)
  ) +
  scale_y_discrete(limits = y_levels, expand = c(0,0))

# ---------- Right strip ----------
label_df <- counts_input %>% distinct(cellType, cellType_short) %>%
  mutate(y = factor(cellType, levels = y_levels))

p_right <- ggplot(label_df, aes(x = 0, y = y, label = cellType_short)) +
  geom_text(hjust = 0, size = 3.0) +
  xlim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(10, 14, 10, 2))

# ---------- Title + Assemble ----------
title_grob <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Pathway Enrichment by Cell Type",
           fontface = "bold", size = 5, hjust = 0) +
  annotate("text", x = 0, y = -0.4,
           label = paste0("Left: counts  |  Right: proportions (sum = 100% per cell type; sig_thr = ", sig_thr, ")"),
           size = 3.4, hjust = 0, colour = "grey30") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 5, 0))

full_plot <- title_grob / (
  p_left + p_counts + p_props + p_right +
    plot_layout(widths = c(1.25, 4.1, 3.7, 1.5), guides = "collect")
)

full_plot

outfile <- file.path(PLOT_DIR, "triptych_counts_props_horizontal.pdf")
ggsave(outfile, full_plot, width = 15.2, height = 9)
message("Saved: ", outfile)


### Stats
counts_input %>%
  group_by(group) %>%
  summarise(total_n = sum(n)) %>%
  mutate(proportion = total_n / sum(total_n))

################################################################################
########################## CONCORDANCE #########################################
################################################################################
############
# Settings #
sig_thr <- 2 # -log10(0.01)
topN <- 10
ct <- "Ast CHI3L1"
#############

concordance_df <- all_fg_wide %>%
  dplyr::transmute(
    pathway,
    cellType,
    mlogq_case = minusLogPAdj_CaseControl,
    mlogq_cpp  = minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD,
    NES_case   = NES_CaseControl,
    NES_cpp    = NES_PathAD.TxAD_vs_PathNonAD.TxNonAD
  ) %>%
  dplyr::mutate(
    score_case = NES_case * mlogq_case,
    score_cpp  = NES_cpp  * mlogq_cpp,
    sig_case = is.finite(mlogq_case) & (mlogq_case >= sig_thr),
    sig_cpp  = is.finite(mlogq_cpp)  & (mlogq_cpp  >= sig_thr),
    group = dplyr::case_when(
      sig_case & sig_cpp ~ "Shared",
      sig_case           ~ "C/C only",
      sig_cpp            ~ "CPP only",
      TRUE               ~ "Not sig"
    ),
    group = factor(group, levels = c("CPP only", "Shared", "C/C only", "Not sig"))
  ) %>%
  dplyr::filter(is.finite(score_case), is.finite(score_cpp))

plot_pathway_concordance <- function(concordance_df, ct, topN = 5) {
  group_colors <- c(
    "CPP only" = "#E41A1C",
    "Shared"   = "#984EA3",
    "C/C only" = "#4DAF4A",
    "Not sig"  = "grey85"
  )
  group_colors <- c(
    "CPP only" = "#e9d985",
    "Shared"   = "#b2bd7e",
    "C/C only" = "#749c75",
    "Not sig"  = "grey85"
  )
  plot_df <- subset(concordance_df, cellType == ct)
  # CPP only
  top_cpp_up <- plot_df %>%
    filter(group == "CPP only", score_cpp > 0) %>%
    arrange(desc(score_cpp)) %>%
    slice_head(n = topN)
  top_cpp_down <- plot_df %>%
    filter(group == "CPP only", score_cpp < 0) %>%
    arrange(score_cpp) %>%
    slice_head(n = topN)
  # C/C only
  top_cc_up <- plot_df %>%
    filter(group == "C/C only", score_case > 0) %>%
    arrange(desc(score_case)) %>%
    slice_head(n = topN)
  top_cc_down <- plot_df %>%
    filter(group == "C/C only", score_case < 0) %>%
    arrange(score_case) %>%
    slice_head(n = topN)
  # Shared (by average score)
  top_shared_up <- plot_df %>%
    filter(group == "Shared") %>%
    mutate(avg_score = (score_case + score_cpp) / 2) %>%
    filter(avg_score > 0) %>%
    arrange(desc(avg_score)) %>%
    slice_head(n = topN)
  top_shared_down <- plot_df %>%
    filter(group == "Shared") %>%
    mutate(avg_score = (score_case + score_cpp) / 2) %>%
    filter(avg_score < 0) %>%
    arrange(avg_score) %>%
    slice_head(n = topN)
  # Combine
  top_paths <- bind_rows(
    top_cpp_up, top_cpp_down,
    top_cc_up, top_cc_down,
    top_shared_up, top_shared_down
  ) %>%
    mutate(pathway_clean = gsub("^(GOBP_|HALLMARK_)", "", pathway),
           pathway_clean = gsub("_", " ", pathway_clean),
           pathway_clean = tools::toTitleCase(tolower(pathway_clean)))
  ggplot(plot_df, aes(x = score_case, y = score_cpp)) +
    # RASTERIZED POINTS - keeps file size small
    ggrastr::geom_point_rast(aes(color = group), 
                             alpha = 0.5, size = 1.5, raster.dpi = 300) +
    ggrastr::geom_point_rast(data = filter(plot_df, group != "Not sig"),
                             aes(color = group), size = 2.5, alpha = 0.8, raster.dpi = 300) +
    # Vector reference lines
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.3) +
    # VECTOR TEXT - stays editable!
    geom_text_repel(data = top_paths, 
                    aes(label = pathway_clean, color = group), 
                    size = 2.8, 
                    max.overlaps = 30,
                    box.padding = 0.5, 
                    point.padding = 0.3,
                    segment.size = 0.2,
                    segment.alpha = 0.6,
                    min.segment.length = 0.1,
                    force = 2,
                    fontface = "bold") +
    scale_color_manual(values = group_colors) +
    labs(
      title = ct,
      x = "NES × -log(q) [Case/Control]",
      y = "NES × -log(q) [CPP]",
      color = "Significance"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      legend.position = "none"
    )
}

topN = 5

excCbln <- plot_pathway_concordance(concordance_df, ct = "Exc L2-3 CBLN2 LINC02306", topN = topN)
excRorb <- plot_pathway_concordance(concordance_df, ct = "Exc L4-5 RORB GABRG1", topN = topN)
excReln <- plot_pathway_concordance(concordance_df, ct = "Exc RELN CHD7", topN = topN)
excThemis <- plot_pathway_concordance(concordance_df, ct = "Exc L6 THEMIS NFIA", topN = topN)
inhVip <- plot_pathway_concordance(concordance_df, ct = "Inh VIP ABI3BP", topN = topN)
astCH <- plot_pathway_concordance(concordance_df, ct = "Ast CHI3L1", topN = topN)
oli <- plot_pathway_concordance(concordance_df, ct = "Oli", topN = topN)
opc <- plot_pathway_concordance(concordance_df, ct = "OPC", topN = topN)
mic <- plot_pathway_concordance(concordance_df, ct = "Mic P2RY12", topN = topN)

pdf("pdf/gsea.Scatter.pdf", height=9, width=24)
grid.arrange(excCbln, excRorb, excReln, excThemis, inhVip, astCH, oli, opc, mic, nrow=2)
dev.off()







################################################################################
################################################################################
################################################################################
############### OLD UNUSED CODE GRAVEYARD BELOW ################################
################################################################################
################################################################################
################################################################################



################################################################################
########################## BIG HEATMAP #########################################
################################################################################

minusLogPAdjThreshold <- 8
fg_heat <- all_fg_wide %>%
  select(pathway, cellType, minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD) %>%
  filter(is.finite(minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD))
fg_heat$minusLogPAdj <- fg_heat$minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD
fg_heat$minusLogPAdj_PathAD.TxAD_vs_PathNonAD.TxNonAD <- NULL

# ---- Filter to pathways significant in at least one cell type ----
sig_paths <- fg_heat %>%
  group_by(pathway) %>%
  filter(any(minusLogPAdj >= minusLogPAdjThreshold)) %>%
  ungroup()

# ---- Pivot to wide format: pathway × cellType matrix ----
heat_mat <- sig_paths %>%
  pivot_wider(
    names_from = cellType,
    values_from = minusLogPAdj,
    values_fill = 0
  )

# Convert to matrix and set row names
mat <- as.matrix(heat_mat[,-1])
rownames(mat) <- heat_mat$pathway

# ---- Heatmap ----
pheatmap(
  mat,
  scale = "none",
  color = colorRampPalette(c("white", "#c1121f"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  main = "Pathway Enrichment (-log10 pAdj) across Cell Types (CPP Contrast)"
)




