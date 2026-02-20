library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(stringr)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(ggpubr)

################################################################################
######################## Making combined file ##################################
################################################################################

whitney_files <- list.files(
  path = "whitney_test",
  pattern = "whitney_fdr05\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

parse_whitney_metadata <- function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  
  tibble(
    cell_type = parts[which(parts == "whitney_test") + 1],
    contrast  = parts[which(parts == "whitney_test") + 2]
  )
}


all_whitney_results <- map_dfr(
  whitney_files,
  function(f) {
    df <- read.csv(f)
    if (nrow(df) == 0) {
      return(NULL)
    }
    meta <- parse_whitney_metadata(f)
    df %>%
      filter(!is.na(q_value)) %>%  # Remove degenerate results
      mutate(
        cell_type = meta$cell_type,
        contrast  = meta$contrast
      )
  }
)

all_whitney_results <- all_whitney_results %>%
  group_by(reaction_name) %>%
  fill(subsystem, compartments_str, reaction_type, metabolites, 
       .direction = "downup") %>%
  ungroup()

subsystemsToExclude <- c(
  "Artificial reactions",
  "Pool reactions",
  "Exchange/demand reactions",
  "Isolated",
  "Miscellaneous",
  "Transport reactions"
)

all_whitney_results <- subset(all_whitney_results, !(subsystem %in% subsystemsToExclude))
write_tsv(all_whitney_results, "all_reaction_results.txt")

all_whitney_results_cpp <- subset(all_whitney_results, contrast == "cpp_based_tx_group")
all_whitney_results_cc <- subset(all_whitney_results, contrast == "path_ad_path_non_ad_group")
all_whitney_results_shuf <- subset(all_whitney_results, contrast == "shuffle_group")

# it seems we don't quite have the same number of cell types across contrasts, we'll progress 
# analysis with the intersection
cell_types_whitney_all <- intersect(
  intersect(unique(all_whitney_results_cpp$cell_type), 
            unique(all_whitney_results_cc$cell_type)),
  unique(all_whitney_results_shuf$cell_type))

all_whitney_results_cpp_filt <- subset(all_whitney_results_cpp, cell_type %in% cell_types_whitney_all)
all_whitney_results_cc_filt <- subset(all_whitney_results_cc, cell_type %in% cell_types_whitney_all)
all_whitney_results_shuf_filt <- subset(all_whitney_results_shuf, cell_type %in% cell_types_whitney_all)

sig_threshold <- 0.05

dim(subset(all_whitney_results_cpp_filt, q_value < sig_threshold))
dim(subset(all_whitney_results_cc_filt, q_value < sig_threshold))
dim(subset(all_whitney_results_shuf_filt, q_value < sig_threshold))

################################################################################
################################ Plotting Counts ###############################
################################################################################
cellOrder <- c(
  "Exc_L3_5_RORB_PLCH1",
  "Exc_L2_3_CBLN2_LINC02306",
  "Exc_L6_THEMIS_NFIA",
  "Exc_L5_6_IT_Car3",
  "Exc_L4_5_RORB_GABRG1",
  "Exc_L3_4_RORB_CUX2",
  "Exc_L5_6_RORB_LINC02196",
  "Exc_L4_5_RORB_IL1RAPL2",
  "Exc_L6_CT",
  "Exc_L6b",
  "Exc_RELN_CHD7",
  "Exc_L5_6_NP",
  "Exc_NRGN",
  "Exc_L5_ET",
  "Inh_L1_6_LAMP5_CA13",
  "Inh_RYR3_TSHZ2",
  "Inh_ALCAM_TRPM3",
  "Inh_SGCD_PDE3A",
  "Inh_LAMP5_NRG1_(Rosehip)",
  "Inh_SORCS1_TTN",  
  "Inh_L6_SST_NPY",
  "Inh_PVALB_SULF1",
  "Inh_L5_6_PVALB_STON2",
  "Inh_VIP_ABI3BP",
  "Inh_PVALB_CA8_(Chandelier)",
  "Inh_L1_2_PAX6_SCGN",
  "Inh_PVALB_HTR4",  
  "Inh_FBN2_EPB41L4A",
  "Inh_PTPRK_FAM19A1",
  "Inh_VIP_CLSTN2",
  "Inh_L5_6_SST_TH",
  "Inh_GPC5_RIT2",
  "Inh_VIP_TSHZ2",
  "Inh_VIP_THSD7B",
  "Inh_ENOX2_SPHKAP",
  "Inh_L3_5_SST_MAFB",
  "Inh_LAMP5_RELN",
  "Inh_L1_PAX6_CA4",
  "Inh_CUX2_MSR1",
  "Ast_GRM3",
  "Ast_DPP10",
  "Ast_CHI3L1",
  "Oli",
  "OPC",
  "Mic_P2RY12",
  "Mic_MKI67",
  "Mic_TPT1",
  "T_cells",
  "CAMs",
  "Per",
  "SMC",
  "End",
  "Fib_FLRT2",
  "Fib_SLC4A4"
)

totalCellTypeCounts <- all_whitney_results_cpp %>%
  count(cell_type, name = "total_reactions_per_cell_type")

# Calculate number of significant reactions per cell type
significantCellTypeCounts <- all_whitney_results_cpp %>%
  filter(q_value < sig_threshold) %>%
  count(cell_type, name = "n_significant_reactions_per_cell_type")

# Merge the two data frames and calculate proportion
cellTypeProportions <- left_join(totalCellTypeCounts, significantCellTypeCounts, by = "cell_type") %>%
  mutate(
    n_significant_reactions_per_cell_type = replace_na(n_significant_reactions_per_cell_type, 0),
    proportion_significant = n_significant_reactions_per_cell_type / total_reactions_per_cell_type
  )

cellTypeProportions$cell_type <- factor(cellTypeProportions$cell_type, levels=rev(cellOrder))

a <- ggplot(cellTypeProportions,
            aes(x = cell_type, y = proportion_significant, fill = cell_type)) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = unit(0.8, "lines")
  ) +
  labs(
    x = "Cell Type",
    y = "Proportion of significant reactions (FDR < 0.05)",
    title = "Proportion of significant reactions per cell type"
  ) 

pdf("pdf/rxnPerCellType.pdf", height=7, width=5)
a
dev.off()

################################################################################
################################ Enrichment ####################################
################################################################################

library(dplyr)
library(tidyr)

enrichment_df <- all_whitney_results_cpp %>%
  group_by(cell_type) %>%
  group_modify(~ {
    df <- .x
    subsystems <- df %>% # find all subsystems significant in at least one cell type
      filter(q_value < sig_threshold) %>%
      pull(subsystem) %>%
      unique()
    
    subsystems %>%
      map_dfr(~ {
        a <- sum(df$subsystem == .x & df$q_value < sig_threshold)
        b <- sum(df$subsystem == .x & df$q_value >= sig_threshold)
        c <- sum(df$subsystem != .x & df$q_value < sig_threshold)
        d <- sum(df$subsystem != .x & df$q_value >= sig_threshold)
        
        mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
        
        if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) {
          tibble(
            subsystem = .x,
            p_enrich = NA_real_,
            odds_ratio = NA_real_,
            a = a, b = b, c = c, d = d
          )
        } else {
          ft <- fisher.test(mat)
          tibble(
            subsystem = .x,
            p_enrich = ft$p.value,
            odds_ratio = unname(ft$estimate),
            a = a, b = b, c = c, d = d
          )
        }
      })
  }) %>%
  ungroup() %>%
  group_by(cell_type) %>%
  mutate(q_enrich = p.adjust(p_enrich, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    minusLogQ = -log10(q_enrich),
    log2OR = log2(odds_ratio)
  )

# Define a custom palette for significant cell types
cell_type_palette_old <- c(
  "Ast_DPP10" = "#e02428",
  "Mic_P2RY12" = "#8a3f9a",
  "Exc_L4_5_RORB_GABRG1" = "#9ef01a",
  "Exc_L2_3_CBLN2_LINC02306" = "#70e000",
  "Exc_L4_5_RORB_IL1RAPL2" = "#38b000",
  "Exc_L5_6_IT_Car3" = "#008000",
  "Exc_L5_ET" = "#006400",
  "Exc_L5_6_RORB_LINC02196" = "#004b23",
  "Inh_PVALB_CA8_(Chandelier)" = "#03045e",
  "Inh_VIP_CLSTN2" = "#023e8a",
  "Inh_L1_6_LAMP5_CA13" = "#0077b6",
  "Inh_L5_6_SST_TH" = "#0096c7",
  "Inh_LAMP5_NRG1_(Rosehip)" = "#00b4d8",
  "Inh_PVALB_HTR4" = "#48cae4",
  "Inh_RYR3_TSHZ2" = "#90e0ef",
  "Inh_VIP_TSHZ2" = "#ade8f4"
)

cell_type_palette <- c(
  "Exc_L3_5_RORB_PLCH1"        = "#9EF01A",
  "Exc_L4_5_RORB_GABRG1"       = "#70E000",
  "Exc_L5_6_IT_Car3"           = "#38B000",
  "Exc_L5_6_RORB_LINC02196"    = "#008000",
  "Exc_L5_ET"                  = "#006400",
  "Exc_L6_THEMIS_NFIA"         = "#004B23",
  "Inh_L5_6_SST_TH"            = "#1F77B4",
  "Inh_PVALB_CA8_(Chandelier)" = "#1B6CA8",
  "Inh_PVALB_HTR4"             = "#155D8B",
  "Inh_RYR3_TSHZ2"             = "#0F4C75",
  "Inh_VIP_CLSTN2"             = "#0A3D62",
  "Inh_VIP_TSHZ2"              = "#082F4E",
  "Mic_P2RY12"                 = "#7B2CBF"
)

# Filter top hits for labeling
top_hits <- enrichment_df %>%
  filter(minusLogQ > 2.3)

# Create the volcano plot
b <- ggplot(enrichment_df,
       aes(x = log2OR,
           y = minusLogQ,
           color = ifelse(minusLogQ > 2, cell_type, "Not Significant"))) +
  geom_point(alpha = 1, size = 3) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = top_hits,
    aes(label = subsystem),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c(cell_type_palette, "Not Significant" = "grey")) +
  theme_classic() +
  labs(
    x = "log2 odds ratio",
    y = "-log10(FDR)",
    title = "Subsystem enrichment of dysregulated reactions",
    color = "Cell Type"
  )

pdf("pdf/subsystemRxnEnrichment.pdf", height=6, width=9)
b
dev.off()
b

write_tsv(enrichment_df, "enrichment_df.txt")
