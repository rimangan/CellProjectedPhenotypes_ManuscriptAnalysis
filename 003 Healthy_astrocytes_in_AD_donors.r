library(Seurat)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(gridExtra)
library(ggridges)

######################################
######## BOXPLOT SUPP FIGURE #########
######################################

cppScores <- readRDS(file = "Data/all_cell_scores_new3.rds")
cppScores <- cppScores %>% mutate(projid = as.integer(as.character(projid)))
Exc <- cppScores$celltype[grepl("^Exc\\.", cppScores$celltype)] |> unique()
Inh <- cppScores$celltype[grepl("^Inh\\.", cppScores$celltype)] |> unique()
Ast <- cppScores$celltype[grepl("^Ast\\.", cppScores$celltype)] |> unique()
Mic <- c("Mic.P2RY12", "Mic.MKI67")
ImmuneVasc <- c("T.cells", "CAMs", "Per", "SMC", "End", "Fib.FLRT2", "Fib.SLC4A4")
Oli <- c("Oli")
OPC <- c("OPCs")
majorCellTypes <- list(
  Exc = Exc,
  Inh = Inh,
  Ast = Ast,
  Mic = Mic,
  ImmuneVasc = ImmuneVasc,
  Oli = Oli,
  OPC = OPC
)

meta <- read.csv("Data/Metadata.csv")
meta$AD_status_donor <- ifelse(meta$niareagansc %in% c(4,3),0,1) ## Binarize AD status
meta <- meta %>% mutate(projid = as.integer(projid))
cppScores <- cppScores %>% left_join(meta %>% select(projid, AD_status_donor), by = "projid")
#cppScores$AD_status_donor <- factor(cppScores$AD_status_donor, levels = c(0, 1))

plot_cpp_boxplot <- function(cppScores, major_celltype, phenotype_col = "AD_status", threshold = 0.25) {
  df <- cppScores %>% filter(celltype %in% major_celltype)
  median_scores <- tapply(df[[phenotype_col]], df$projid, median, na.rm=TRUE)
  df$projid <- factor(df$projid, levels = names(sort(median_scores)))
  # Plot
  p <- ggplot(df, aes(x = projid, y = .data[[phenotype_col]], fill = AD_status_donor, color=AD_status_donor)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = c(-threshold, threshold), linetype = "dashed") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none"
          ) +
    scale_fill_gradientn(colors = c("blue", "red"), limits = c(0, 1)) +
    scale_color_gradientn(colors = c("blue", "red"), limits = c(0, 1)) +
    labs(
      y = "CPP Score"
    )
  return(p)
}


make_column <- function(cppScores, currPhenotype) {
  pExc <- plot_cpp_boxplot(cppScores, Exc, currPhenotype) + labs(title = "Excitatory Neurons")
  pInh <- plot_cpp_boxplot(cppScores, Inh, currPhenotype) + labs(title = "Inhibitory Neurons")
  pAst <- plot_cpp_boxplot(cppScores, Ast, currPhenotype) + labs(title = "Astrocytes")
  pMic <- plot_cpp_boxplot(cppScores, Mic, currPhenotype) + labs(title = "Microglia")
  pImmuneVasc <- plot_cpp_boxplot(cppScores, ImmuneVasc, currPhenotype) + labs(title = "Immune & Vascular Cells")
  pOli <- plot_cpp_boxplot(cppScores, Oli, currPhenotype) + labs(title = "Oligodendrocytes")
  pOPC <- plot_cpp_boxplot(cppScores, OPC, currPhenotype) + labs(title = "OPCs")
  return(pExc / pInh / pAst / pMic / pImmuneVasc / pOli / pOPC)
}


row_AD <- make_column(cppScores, "AD_status")
row_amyloid <- make_column(cppScores, "amyloid")
row_tangles <- make_column(cppScores, "tangles")

final_plot <- row_AD | row_amyloid | row_tangles
pdf("Figures/Supplementary/Supp3_Boxplots/allAdPlots.pdf", height=10, width=12)
final_plot
dev.off()

pdf("Figures/Supplementary/Supp3_Boxplots/astAD_status_Zoom.pdf", height=4, width=8)
plot_cpp_boxplot(cppScores, Ast, "AD_status") + labs(x="377 donors sorted by median cell projected phenotype score for AD status", y="CPP Score (AD Status)")
dev.off()

#####################################
### Figure S3a - but as a ridge #####
#####################################
plot_cpp_ridge <- function(cppScores, major_celltype, phenotype_col = "AD_status", 
                          threshold = 0.25, max_donors = 100, seed = 123) {
  # Filter for major cell types
  df <- cppScores %>% filter(celltype %in% major_celltype)
  
  # Randomly downsample donors if needed
  set.seed(seed)
  donors <- unique(df$projid)
  if(length(donors) > max_donors) {
    sampled_donors <- sample(donors, max_donors)
    df <- df %>% filter(projid %in% sampled_donors)
  }
  
  # Order donors by median score
  median_scores <- tapply(df[[phenotype_col]], df$projid, median, na.rm = TRUE)
  df$projid <- factor(df$projid, levels = names(sort(median_scores)))
  
  # Ridge plot
  p <- ggplot(df, aes(x = .data[[phenotype_col]], y = projid, fill = AD_status_donor, color = AD_status_donor)) +
    geom_density_ridges(scale = 1, rel_min_height = 0.01, alpha = 0.8) +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_gradientn(colors = c("blue", "red"), limits = c(0, 1)) +
    scale_color_gradientn(colors = c("blue", "red"), limits = c(0, 1)) +
    labs(
      x = "CPP Score"
    )
  
  return(p)
}

p <- plot_cpp_ridge(cppScores, Mic, max_donors = 50)
p + labs(x="AD Status CPP Score (Microglia)", y="50 Postmortem Donors") + theme(axis.text.y=element_blank(), axis.title.y=element_text(size=12, angle=90))

#####################################
#### Figure S3b - stacked bar #######
#####################################

cppScoresAst <- filter(cppScores, celltype %in% Ast)
#count AD cells
cppScoresAst %>%
  filter(AD_status > 0.25) %>%
  count()
#count nonAD cells
cppScoresAst %>%
  filter(AD_status < -0.25) %>%
  count()
#count neutral cells
cppScoresAst %>%
  filter(AD_status >-0.25 & AD_status < 0.25) %>%
  count()
# disease from AD donors
cppScoresAst %>%
  filter(AD_status > 0.25 & AD_status_donor == 1) %>%
  count()
#healthy from non AD donors
cppScoresAst %>%
  filter(AD_status < -0.25 & AD_status_donor == 0) %>%
  count()
#healthy from AD donors
cppScoresAst %>%
  filter(AD_status < -0.25 & AD_status_donor == 1) %>%
  count()
#disease from nonAD donors
cppScoresAst %>%
  filter(AD_status > 0.25 & AD_status_donor == 0) %>%
  count()

cppScoresAst <- cppScoresAst %>%
  mutate(cell_state = case_when(
    AD_status > 0.25 ~ "AD-like",
    AD_status < -0.25 ~ "Healthy-like",
    TRUE ~ "Neutral"
  ),
  AD_status_donor = ifelse(AD_status_donor == 1, "AD", "NonAD"))
cppScoresAst$cell_state <- factor(
  cppScoresAst$cell_state,
  levels = c("AD-like", "Neutral", "Healthy-like")
)
# Plot of cell states by donor
df_counts <- cppScoresAst %>%
  count(AD_status_donor, cell_state)
df_counts <- df_counts %>%
  group_by(AD_status_donor) %>%
  mutate(percent = n / sum(n) * 100)
pdf("Figures/Supplementary/Supp3_Boxplots/stateDistribution.pdf", height=6, width=6)
ggplot(df_counts, aes(x = AD_status_donor, y = percent, fill = cell_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("AD-like" = "#ec2426", "Neutral" = "#231f20", "Healthy-like" = "#3b54a4")) +
  labs(
    x = "Donor Clinical Diagnosis",
    y = "Percent of Astrocytes",
    fill = "CPP-Inferred Cell State",
    title = "Distribution of Astrocyte Cell States by Donor AD Status"
  ) +
  theme_classic()
dev.off()


####### for ASHG presentation ###############
cppScoresMic <- filter(cppScores, celltype %in% "Mic.P2RY12")
# Count AD-like cells
cppScoresMic %>%
  filter(AD_status > 0.25) %>%
  count()

# Count non-AD-like cells
cppScoresMic %>%
  filter(AD_status < -0.25) %>%
  count()

# Count neutral cells
cppScoresMic %>%
  filter(AD_status > -0.25 & AD_status < 0.25) %>%
  count()

# AD-like cells from AD donors
cppScoresMic %>%
  filter(AD_status > 0.25 & AD_status_donor == 1) %>%
  count()

# Healthy-like cells from non-AD donors
cppScoresMic %>%
  filter(AD_status < -0.25 & AD_status_donor == 0) %>%
  count()

# Healthy-like cells from AD donors
cppScoresMic %>%
  filter(AD_status < -0.25 & AD_status_donor == 1) %>%
  count()

# AD-like cells from non-AD donors
cppScoresMic %>%
  filter(AD_status > 0.25 & AD_status_donor == 0) %>%
  count()

cppScoresMic <- cppScoresMic %>%
  mutate(cell_state = case_when(
    AD_status > 0.25 ~ "AD-like",
    AD_status < -0.25 ~ "Healthy-like",
    TRUE ~ "Neutral"
  ),
  AD_status_donor = ifelse(AD_status_donor == 1, "AD", "NonAD"))
cppScoresMic$cell_state <- factor(
  cppScoresMic$cell_state,
  levels = c("Healthy-like", "Neutral", "AD-like")
)

df_counts <- cppScoresMic %>%
  count(AD_status_donor, cell_state)
df_counts <- df_counts %>%
  group_by(AD_status_donor) %>%
  mutate(percent = n / sum(n) * 100)
df_counts$AD_status_donor <- factor(
  df_counts$AD_status_donor,
  levels = c("NonAD", "AD")
)

# Plot distribution of microglial states by donor diagnosis
pdf("pdf/ashg_rotated_panel.pdf", height = 4, width = 9)
ggplot(df_counts, aes(x = AD_status_donor, y = percent, fill = cell_state)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "AD-like" = "#ec2426",
    "Neutral" = "#231f20",
    "Healthy-like" = "#3b54a4"
  )) +
  labs(
    x = "Donor Clinical Diagnosis",
    y = "Percent of Microglia",
    fill = "CPP-Inferred Cell State",
    title = "Distribution of Microglial Cell States by Donor AD Status"
  ) +
  theme_classic() +
  coord_flip()
dev.off()

pdf("pdf/mic.barPlot.Ad_status.pdf", height=6, width=9)
plot_cpp_boxplot(cppScores, Mic, "AD_status") + labs(x="377 donors sorted by median cell projected phenotype score for AD status", y="CPP Score (AD Status)")
dev.off()

library(tidyverse)
library(ggridges)

# Pick 2 AD donors and 2 non-AD donors at random
set.seed(42)  # for reproducibility

donor_subset <- cppScores %>%
  filter(celltype %in% Mic) %>%
  group_by(AD_status_donor) %>%
  sample_n(2) %>%
  pull(projid) %>%
  unique()

donor_subset

# Subset data for ridge plot
df_subset <- cppScores %>%
  filter(celltype %in% Mic, projid %in% donor_subset) %>%
  mutate(
    projid = factor(projid, levels = unique(projid)),  # ensure categorical
    AD_status_donor = factor(AD_status_donor, levels = c(0, 1), labels = c("Non-AD", "AD"))
  )

# Ridge plot: distributions of CPP scores for the selected donors
p_zoom <- ggplot(df_subset, aes(x = AD_status, y = projid, fill = as.factor(AD_status_donor))) +
  geom_density_ridges(alpha = 0.7, scale = 1.1, rel_min_height = 0.01, color = "white") +
  scale_fill_manual(
    values = c("0" = "#3b54a4", "1" = "#ec2426"),
    labels = c("Non-AD donor", "AD donor")
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray40") +
  labs(
    x = "CPP Score (AD Status)",
    y = "Donor ID",
    fill = "Donor Diagnosis",
    title = "Zoomed-in distributions of microglial CPP scores across selected donors"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 10)
  )

p_zoom



#####################################
###### DEG Analysis on Mic ##########
#####################################

#if running just this bit in a fresh session, run these lines:
cppScores <- readRDS(file = "Data/all_cell_scores_new3.rds")
if (!"cellId" %in% colnames(cppScores)) {
  cppScores$cellId <- rownames(cppScores)
}
cppScores <- cppScores %>% mutate(projid = as.integer(as.character(projid)))
Mic <- cppScores$celltype[grepl("^Mic\\.", cppScores$celltype)] |> unique()
meta <- read.csv("Data/Metadata.csv")
meta$AD_status_donor <- ifelse(meta$niareagansc %in% c(4,3),0,1) ## Binarize AD status
meta <- meta %>% mutate(projid = as.integer(projid))
cppScores <- cppScores %>% left_join(meta %>% select(projid, AD_status_donor), by = "projid")
cppScoresMic <- filter(cppScores, celltype %in% Mic)
rownames(cppScoresMic) <- cppScoresMic$cellId

#otherwise, start here
mic_cells <- readRDS("~/KellisProjects/CPP/GerardAD430/CPP/Data/mathys/Immune_cells.rds")
mic_cells <- UpdateSeuratObject(mic_cells) # if running on seurat 5. This rds was made in seurat 4
mic_cells <- AddMetaData(mic_cells, metadata = cppScoresMic)
mic_cells <- NormalizeData(mic_cells)

# here we compare healthy cells in Ad vs. non-AD donors
healthy_cells <- WhichCells(mic_cells, expression = AD_status < -0.25)
astro_cells_healthy <- subset(mic_cells, cells = healthy_cells)
## Perform DEA between AD and non-AD, but only using their "healthy-like" cells
adVsNonAd_healthy_cells <- FindMarkers(astro_cells_healthy, group.by = "AD_status_donor",ident.1 = "0",logfc.threshold = 0.1)
table(adVsNonAd_healthy_cells$p_val_adj < 0.01)


## here i compare healthy cells from non-AD donors to disease cells from AD donors
disease_cells_from_AD <- WhichCells(mic_cells, expression = AD_status > 0.25 & AD_status_donor == "1")
healthy_nonad_cells <- WhichCells(mic_cells, expression = AD_status < -0.25 & AD_status_donor == "0")
compare_cells <- c(disease_cells_from_AD, healthy_nonad_cells)
mic_cells_subset <- subset(mic_cells, cells = compare_cells)
mic_cells_subset$compare_group <- ifelse(
  mic_cells_subset$AD_status > 0.25, "disease_AD", "healthy_nonAD"
)
diseaseVsHealthy <- FindMarkers(
  mic_cells_subset,
  group.by = "compare_group",
  ident.1 = "disease_AD",
  ident.2 = "healthy_nonAD",
  logfc.threshold = 0.1
)

table(adVsNonAd_healthy_cells$p_val_adj < 0.01)
table(diseaseVsHealthy$p_val_adj < 0.01)

1494/40

write.csv(adVsNonAd_healthy_cells, "Results/SuppTables/STable1.adVsNonAd_healthy_cells.csv")
write.csv(diseaseVsHealthy, "Results/SuppTables/STable2.diseaseVsHealthy.csv")


#AD-like cells from AD donors vs AD-like cells from non-AD donors
adlike_from_AD <- WhichCells(mic_cells, expression = AD_status > 0.25 & AD_status_donor == "1")
adlike_from_nonAD <- WhichCells(mic_cells, expression = AD_status > 0.25 & AD_status_donor == "0")
compare_cells_1 <- c(adlike_from_AD, adlike_from_nonAD)
mic_cells_subset_1 <- subset(mic_cells, cells = compare_cells_1)
mic_cells_subset_1$compare_group <- ifelse(
  mic_cells_subset_1$AD_status_donor == "1", "adlike_AD", "adlike_nonAD"
)
adlike_ADvsNonAD <- FindMarkers(
  mic_cells_subset_1,
  group.by = "compare_group",
  ident.1 = "adlike_AD",
  ident.2 = "adlike_nonAD",
  logfc.threshold = 0.1
)

# non-AD-like cells from AD donors vs AD-like cells from non-AD donors
nonadlike_from_AD <- WhichCells(mic_cells, expression = AD_status < -0.25 & AD_status_donor == "1")
adlike_from_nonAD <- WhichCells(mic_cells, expression = AD_status > 0.25 & AD_status_donor == "0")
compare_cells_2 <- c(nonadlike_from_AD, adlike_from_nonAD)
mic_cells_subset_2 <- subset(mic_cells, cells = compare_cells_2)
mic_cells_subset_2$compare_group <- ifelse(
  mic_cells_subset_2$AD_status < -0.25, "nonadlike_AD", "adlike_nonAD"
)
nonadlikeAD_vs_adlikeNonAD <- FindMarkers(
  mic_cells_subset_2,
  group.by = "compare_group",
  ident.1 = "nonadlike_AD",
  ident.2 = "adlike_nonAD",
  logfc.threshold = 0.1
)

table(adlike_ADvsNonAD$p_val_adj < 0.01)
table(nonadlikeAD_vs_adlikeNonAD$p_val_adj < 0.01)

#AD-like vs non-AD-like cells, within AD donors only
adlike_from_AD <- WhichCells(mic_cells, expression = AD_status > 0.25 & AD_status_donor == "1")
nonadlike_from_AD <- WhichCells(mic_cells, expression = AD_status < -0.25 & AD_status_donor == "1")
compare_cells_3 <- c(adlike_from_AD, nonadlike_from_AD)
mic_cells_subset_3 <- subset(mic_cells, cells = compare_cells_3)
mic_cells_subset_3$compare_group <- ifelse(
  mic_cells_subset_3$AD_status > 0.25, "adlike_AD", "nonadlike_AD"
)
withinAD_adlike_vs_nonadlike <- FindMarkers(
  mic_cells_subset_3,
  group.by = "compare_group",
  ident.1 = "adlike_AD",
  ident.2 = "nonadlike_AD",
  logfc.threshold = 0.1
)

# AD-like vs non-AD-like cells, within non-AD donors only
adlike_from_nonAD <- WhichCells(mic_cells, expression = AD_status > 0.25 & AD_status_donor == "0")
nonadlike_from_nonAD <- WhichCells(mic_cells, expression = AD_status < -0.25 & AD_status_donor == "0")
compare_cells_4 <- c(adlike_from_nonAD, nonadlike_from_nonAD)
mic_cells_subset_4 <- subset(mic_cells, cells = compare_cells_4)
mic_cells_subset_4$compare_group <- ifelse(
  mic_cells_subset_4$AD_status > 0.25, "adlike_nonAD", "nonadlike_nonAD"
)
withinNonAD_adlike_vs_nonadlike <- FindMarkers(
  mic_cells_subset_4,
  group.by = "compare_group",
  ident.1 = "adlike_nonAD",
  ident.2 = "nonadlike_nonAD",
  logfc.threshold = 0.1
)

table(diseaseVsHealthy$p_val_adj < 0.01)             # ADdonor:ADlike vs NonADdonor:NonADlike (diagonal) table s1
table(adVsNonAd_healthy_cells$p_val_adj < 0.01)      # ADdonor:NonADlike vs NonADdonor:NonADlike table s6
table(adlike_ADvsNonAD$p_val_adj < 0.01)              # ADdonor:ADlike vs NonADdonor:ADlike table s5
table(nonadlikeAD_vs_adlikeNonAD$p_val_adj < 0.01)   # ADdonor:NonADlike vs NonADdonor:ADlike (anti-diagonal) table s4
table(withinAD_adlike_vs_nonadlike$p_val_adj < 0.01)  # ADdonor:ADlike vs ADdonor:NonADlike table s2
table(withinNonAD_adlike_vs_nonadlike$p_val_adj < 0.01) # NonADdonor:ADlike vs NonADdonor:NonADlike table s3

# Table S1: AD-like cells from AD donors vs Control-like cells from non-AD donors (concordant extremes)
write.csv(diseaseVsHealthy, "Tables/TableS1.ADlike_ADdonor_vs_Controllike_NonADdonor.csv")

# Table S2: AD-like vs Control-like cells within AD donors
write.csv(withinAD_adlike_vs_nonadlike, "Tables/TableS2.ADlike_vs_Controllike_withinAD.csv")

# Table S3: AD-like vs Control-like cells within non-AD donors
write.csv(withinNonAD_adlike_vs_nonadlike, "Tables/TableS3.ADlike_vs_Controllike_withinNonAD.csv")

# Table S4: Control-like cells from AD donors vs AD-like cells from non-AD donors (discordant)
write.csv(nonadlikeAD_vs_adlikeNonAD, "Tables/TableS4.Controllike_ADdonor_vs_ADlike_NonADdonor.csv")

# Table S5: AD-like cells from AD donors vs AD-like cells from non-AD donors
write.csv(adlike_ADvsNonAD, "Tables/TableS5.ADlike_ADdonor_vs_ADlike_NonADdonor.csv")

# Table S6: Control-like cells from AD donors vs Control-like cells from non-AD donors
write.csv(adVsNonAd_healthy_cells, "Tables/TableS6.Controllike_ADdonor_vs_Controllike_NonADdonor.csv")

################################################################################
######################## Jaccard Plots #########################################
# Define desired order
contrast_order <- c(
  "ADl_ADd_vs_NonADl_NonADd",    # 1. AD-like:ADdonors vs Control-like:Non-ADdonors
  "ADl_ADd_vs_NonADl_ADd",       # 2. AD-like:ADdonors vs Control-like:ADdonors
  "ADl_NonADd_vs_NonADl_NonADd", # 3. AD-like:nonADdonors vs Control-like:nonADdonors
  "NonADl_ADd_vs_ADl_NonADd",    # 4. Control-like:ADdonors vs AD-like:nonADdonors
  "ADl_ADd_vs_ADl_NonADd",       # 5. AD-like:ADdonors vs AD-like:nonADdonors
  "NonADl_ADd_vs_NonADl_NonADd"  # 6. Control-like:ADdonors vs Control-like:nonADdonors
)

results_list <- list(
  "ADl_ADd_vs_NonADl_NonADd"    = diseaseVsHealthy,
  "ADl_ADd_vs_NonADl_ADd"       = withinAD_adlike_vs_nonadlike,
  "ADl_NonADd_vs_NonADl_NonADd" = withinNonAD_adlike_vs_nonadlike,
  "NonADl_ADd_vs_ADl_NonADd"    = nonadlikeAD_vs_adlikeNonAD,
  "ADl_ADd_vs_ADl_NonADd"       = adlike_ADvsNonAD,
  "NonADl_ADd_vs_NonADl_NonADd" = adVsNonAd_healthy_cells
)

deg_up   <- lapply(results_list, get_sig_up)
deg_down <- lapply(results_list, get_sig_down)
deg_all  <- lapply(results_list, get_sig_all)

jaccard_up   <- compute_jaccard_mat(deg_up)
jaccard_down <- compute_jaccard_mat(deg_down)
jaccard_all  <- compute_jaccard_mat(deg_all)

n <- length(results_list)
nms <- names(results_list)
jaccard_concordant <- matrix(0, n, n, dimnames = list(nms, nms))
jaccard_discordant <- matrix(0, n, n, dimnames = list(nms, nms))

for (i in 1:n) {
  for (j in 1:n) {
    shared_all <- intersect(deg_all[[i]], deg_all[[j]])
    if (length(shared_all) == 0) {
      jaccard_concordant[i, j] <- 0
      jaccard_discordant[i, j] <- 0
      next
    }
    lfc_i <- setNames(results_list[[i]]$avg_log2FC, rownames(results_list[[i]]))
    lfc_j <- setNames(results_list[[j]]$avg_log2FC, rownames(results_list[[j]]))
    same_dir <- shared_all[sign(lfc_i[shared_all]) == sign(lfc_j[shared_all])]
    opp_dir  <- shared_all[sign(lfc_i[shared_all]) != sign(lfc_j[shared_all])]
    all_union <- union(deg_all[[i]], deg_all[[j]])
    jaccard_concordant[i, j] <- length(same_dir) / length(all_union)
    jaccard_discordant[i, j] <- length(opp_dir) / length(all_union)
  }
}

p3 <- pheatmap(jaccard_concordant, display_numbers = TRUE, number_format = "%.2f",
               color = colorRampPalette(c("white", "forestgreen"))(50),
               main = "Jaccard: Concordant (same direction)", fontsize = 7, 
               silent = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)

p4 <- pheatmap(jaccard_discordant, display_numbers = TRUE, number_format = "%.2f",
               color = colorRampPalette(c("white", "darkorange"))(50),
               main = "Jaccard: Discordant (opposite direction)", fontsize = 7, 
               silent = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)


asymmetric_lfc_cor <- function(results_list, threshold = 0.01) {
  n <- length(results_list)
  nms <- names(results_list)
  cor_mat <- matrix(NA, n, n, dimnames = list(nms, nms))
  pval_mat <- matrix(NA, n, n, dimnames = list(nms, nms))
  ngene_mat <- matrix(0, n, n, dimnames = list(nms, nms))
  
  for (i in 1:n) {
    # Significant DEGs from contrast i
    sig_genes_i <- rownames(results_list[[i]])[results_list[[i]]$p_val_adj < threshold]
    lfc_i <- setNames(results_list[[i]]$avg_log2FC, rownames(results_list[[i]]))
    
    for (j in 1:n) {
      lfc_j <- setNames(results_list[[j]]$avg_log2FC, rownames(results_list[[j]]))
      
      # Genes significant in i AND tested in j
      shared <- intersect(sig_genes_i, names(lfc_j))
      ngene_mat[i, j] <- length(shared)
      
      if (length(shared) < 3) {
        cor_mat[i, j] <- NA
        pval_mat[i, j] <- NA
        next
      }
      
      test <- cor.test(lfc_i[shared], lfc_j[shared], method = "spearman", exact = FALSE)
      cor_mat[i, j] <- test$estimate
      pval_mat[i, j] <- test$p.value
    }
  }
  
  list(cor = cor_mat, pval = pval_mat, ngenes = ngene_mat)
}

asym <- asymmetric_lfc_cor(results_list)

display_mat <- matrix(
  #paste0(sprintf("%.2f", asym$cor), pval_to_stars(asym$pval), "\n(n=", asym$ngenes, ")"),
  paste0(sprintf("%.2f", asym$cor)),
  nrow = nrow(asym$cor),
  dimnames = dimnames(asym$cor)
)

c <- pheatmap(asym$cor,
              display_numbers = display_mat,
              number_color = "black",
              fontsize_number = 7,
              color = colorRampPalette(c("steelblue", "white", "firebrick"))(50),
              breaks = seq(-1, 1, length.out = 51),
              main = "Asymmetric log2FC Correlation",
              cluster_rows = FALSE, cluster_cols = FALSE,
              fontsize = 8, silent = TRUE)

pdf("pdf/jaccardAndCor.microgliaMarkers.pdf", height=5, width=15)
grid.arrange(p3[[4]], p4[[4]], c[[4]], nrow = 1)
dev.off()









