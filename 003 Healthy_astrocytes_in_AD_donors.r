library(Seurat)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(gridExtra)
library(ggridges)

######################################
############## SETUP #################
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



### BOXPLOTS ###

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

filterCppScores <- function(cppScores, max_donors = 100, seed = 123) {
  set.seed(seed)
  donors <- unique(cppScores$projid)
  
  if (length(donors) > max_donors) {
    sampled <- sample(donors, max_donors)
    cppScores <- cppScores %>% filter(projid %in% sampled)
  }
  
  return(cppScores)
}

plot_cpp_ridge <- function(cppScores_sampled, major_celltype, 
                           phenotype_col = "AD_status",
                           threshold = 0.25) {
  
  # Filter for major cell types
  df <- cppScores_sampled %>% filter(celltype %in% major_celltype)
  
  # Order donors by median CPP score
  median_scores <- tapply(df[[phenotype_col]], df$projid, median, na.rm = TRUE)
  df$projid <- factor(df$projid, levels = names(sort(median_scores)))
  
  # Ridge plot with heavy overlap
  p <- ggplot(df, aes(
    x = .data[[phenotype_col]], 
    y = projid, 
    fill = AD_status_donor, 
    color = AD_status_donor
  )) +
    geom_density_ridges(
      scale = 3.5,               # make ridges much closer / overlapping
      rel_min_height = 0.01,
      alpha = 0.6,               # transparency for overlapping ridges
      size = 0.3
    ) +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_gradientn(colors = c("blue", "red"), limits = c(0, 1)) +
    scale_color_gradientn(colors = c("blue", "red"), limits = c(0, 1)) +
    labs(x = "CPP Score")
  
  return(p)
}

cppScores_sampled <- filterCppScores(cppScores, max_donors = 100)
pdf("pdf/ridgePlots/mic.AllDonors.ridge.pdf", height=16, width=6)
plot_cpp_ridge(cppScores, major_celltype = Mic)
dev.off()

excAD <- plot_cpp_ridge(cppScores_sampled, major_celltype = Exc) + xlim(-1,1)
inhAD <- plot_cpp_ridge(cppScores_sampled, major_celltype = Inh) + xlim(-1,1)
astAD <- plot_cpp_ridge(cppScores_sampled, major_celltype = Ast) + xlim(-1,1)
immuneVascAD <- plot_cpp_ridge(cppScores_sampled, major_celltype = ImmuneVasc) + xlim(-1,1)
micAD <- plot_cpp_ridge(cppScores_sampled, major_celltype = Mic) + xlim(-1,1)
OliAD <- plot_cpp_ridge(cppScores_sampled, major_celltype = Oli) + xlim(-1,1)
OpcAD <- plot_cpp_ridge(cppScores_sampled, major_celltype = OPC) + xlim(-1,1)

excAmyloid <- plot_cpp_ridge(cppScores_sampled, major_celltype = Exc, phenotype_col = "amyloid") + xlim(-1,1)
inhAmyloid <- plot_cpp_ridge(cppScores_sampled, major_celltype = Inh, phenotype_col = "amyloid") + xlim(-1,1)
astAmyloid <- plot_cpp_ridge(cppScores_sampled, major_celltype = Ast, phenotype_col = "amyloid") + xlim(-1,1)
immuneVascAmyloid <- plot_cpp_ridge(cppScores_sampled, major_celltype = ImmuneVasc, phenotype_col = "amyloid") + xlim(-1,1)
micAmyloid <- plot_cpp_ridge(cppScores_sampled, major_celltype = Mic, phenotype_col = "amyloid") + xlim(-1,1)
oliAmyloid <- plot_cpp_ridge(cppScores_sampled, major_celltype = Oli, phenotype_col = "amyloid") + xlim(-1,1)
opcAmyloid <- plot_cpp_ridge(cppScores_sampled, major_celltype = OPC, phenotype_col = "amyloid") + xlim(-1,1)

excTangles <- plot_cpp_ridge(cppScores_sampled, major_celltype = Exc, phenotype_col = "tangles") + xlim(-1,1)
inhTangles <- plot_cpp_ridge(cppScores_sampled, major_celltype = Inh, phenotype_col = "tangles") + xlim(-1,1)
astTangles <- plot_cpp_ridge(cppScores_sampled, major_celltype = Ast, phenotype_col = "tangles") + xlim(-1,1)
immuneVascTangles <- plot_cpp_ridge(cppScores_sampled, major_celltype = ImmuneVasc, phenotype_col = "tangles") + xlim(-1,1)
micTangles <- plot_cpp_ridge(cppScores_sampled, major_celltype = Mic, phenotype_col = "tangles") + xlim(-1,1)
oliTangles <- plot_cpp_ridge(cppScores_sampled, major_celltype = Oli, phenotype_col = "tangles") + xlim(-1,1)
opcTangles <- plot_cpp_ridge(cppScores_sampled, major_celltype = OPC, phenotype_col = "tangles") + xlim(-1,1)

pdf("pdf/ridgePlots/ridgeGrid.pdf", height=12, width=10)
grid.arrange(
  excAD, inhAD, astAD, immuneVascAD, micAD, OliAD, OpcAD,
  excAmyloid, inhAmyloid, astAmyloid, immuneVascAmyloid, micAmyloid, oliAmyloid, opcAmyloid,
  excTangles, inhTangles, astTangles, immuneVascTangles, micTangles, oliTangles, opcTangles,
  nrow=3
)
dev.off()




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
cppScoresMic <- cppScores %>%
  filter(celltype %in% c("Mic.P2RY12", "Mic.MKI67", "Mic.TPT1")) %>%
  mutate(
    cell_state = case_when(
      AD_status >  0.25 ~ "AD-like",
      AD_status < -0.25 ~ "Healthy-like",
      TRUE              ~ "Neutral"
    ),
    cell_state = factor(cell_state, levels = c("Healthy-like", "Neutral", "AD-like"))
  )

df_counts <- cppScoresMic %>%
  count(AD_status_donor, cell_state) %>%
  group_by(AD_status_donor) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup()

# Plot distribution of microglial states by donor diagnosis
pdf("pdf/ridgePlots/rotated_proportion_panel.pdf", height = 4, width = 9)
ggplot(df_counts, aes(
  x = factor(AD_status_donor, labels = c("Control", "AD")),
  y = percent,
  fill = cell_state
)) +
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
###### DEG Analysis on Ast ##########
#####################################

#if running just this bit in a fresh session, run these lines:
cppScores <- readRDS(file = "Data/all_cell_scores_new3.rds")
if (!"cellId" %in% colnames(cppScores)) {
  cppScores$cellId <- rownames(cppScores)
}
cppScores <- cppScores %>% mutate(projid = as.integer(as.character(projid)))
Ast <- cppScores$celltype[grepl("^Ast\\.", cppScores$celltype)] |> unique()
meta <- read.csv("Data/Metadata.csv")
meta$AD_status_donor <- ifelse(meta$niareagansc %in% c(4,3),0,1) ## Binarize AD status
meta <- meta %>% mutate(projid = as.integer(projid))
cppScores <- cppScores %>% left_join(meta %>% select(projid, AD_status_donor), by = "projid")
cppScoresAst <- filter(cppScores, celltype %in% Ast)
rownames(cppScoresAst) <- cppScoresAst$cellId

#otherwise, start here
astro_cells <- readRDS("Data/mathys/Astrocytes.rds")
astro_cells <- UpdateSeuratObject(astro_cells) # if running on seurat 5. This rds was made in seurat 4
astro_cells <- AddMetaData(astro_cells, metadata = cppScoresAst)
astro_cells <- NormalizeData(astro_cells)

# here we compare healthy cells in Ad vs. non-AD donors
healthy_cells <- WhichCells(astro_cells, expression = AD_status < -0.25)
astro_cells_healthy <- subset(astro_cells, cells = healthy_cells)
## Perform DEA between AD and non-AD, but only using their "healthy-like" cells
adVsNonAd_healthy_cells <- FindMarkers(astro_cells_healthy, group.by = "AD_status_donor",ident.1 = "0",logfc.threshold = 0.1)
table(adVsNonAd_healthy_cells$p_val_adj < 0.01)


## here i compare healthy cells from non-AD donors to disease cells from AD donors
disease_cells_from_AD <- WhichCells(astro_cells, expression = AD_status > 0.25 & AD_status_donor == "1")
healthy_nonad_cells <- WhichCells(astro_cells, expression = AD_status < -0.25 & AD_status_donor == "0")
compare_cells <- c(disease_cells_from_AD, healthy_nonad_cells)
astro_cells_subset <- subset(astro_cells, cells = compare_cells)
astro_cells_subset$compare_group <- ifelse(
  astro_cells_subset$AD_status > 0.25, "disease_AD", "healthy_nonAD"
)
diseaseVsHealthy <- FindMarkers(
  astro_cells_subset,
  group.by = "compare_group",
  ident.1 = "disease_AD",
  ident.2 = "healthy_nonAD",
  logfc.threshold = 0.1
)

table(adVsNonAd_healthy_cells$p_val_adj < 0.01)
table(diseaseVsHealthy$p_val_adj < 0.01)

3571/380

write.csv(adVsNonAd_healthy_cells, "Results/SuppTables/STable1.adVsNonAd_healthy_cells.csv")
write.csv(diseaseVsHealthy, "Results/SuppTables/STable2.diseaseVsHealthy.csv")


















##### OLD (from Gerard) ######
library(Seurat)
library(ggplot2)
library(openxlsx)

#Load meta data
meta <- read.csv("Data/Metadata.csv")
## Binarize AD status
meta$AD_status <- ifelse(meta$niareagansc %in% c(4,3),0,1)
## Calculate Cognitive resilience 
meta$CR <- lm(cogn_global_lv~gpath,data = meta)$residuals

## Load original scRNAseq data
raw_files <- c("Data/mathys/Astrocytes.rds")
raw_cells <- readRDS(sprintf("./%s",raw_files))

## Load cell scores
allcombined <- readRDS(file = "Data/all_cell_scores_new3.rds")


## Load cell meta data
astroCell_metadata <- lapply(raw_files,function(file){
    message(file)
    raw_cells <- readRDS(sprintf("./%s",file))
    out <- cbind(raw_cells@reductions$umap@cell.embeddings[,c(1:2)],
                 as.character(raw_cells@meta.data$cell_type_high_resolution)) |> as.data.frame()
    colnames(out) <- c("UMAP_1","UMAP_2","celltype")
    out$main <- gsub(".rds","",file)
    return(out)
})
astroCell_metadata <- astroCell_metadata[[1]]
astroCell_metadata$UMAP_1 <- as.numeric(astroCell_metadata$UMAP_1)
astroCell_metadata$UMAP_2 <- as.numeric(astroCell_metadata$UMAP_2)

astroCell_metadata <- astroCell_metadata[rownames(astroCell_metadata) %in% rownames(allcombined),]
allcombined <- allcombined[rownames(allcombined) %in% rownames(astroCell_metadata),]

astroCell_metadata <- astroCell_metadata[rownames(allcombined),]
astroCell_metadata$projid <- allcombined$projid
astroCell_metadata$AD_status <- allcombined$AD_status

astroCell_metadata$AD_status_individual <- meta[match(astroCell_metadata$projid,meta$projid),"AD_status"] |> as.factor()
astroCell_metadata$cellstatus <- ifelse(astroCell_metadata$AD_status<=-0.25,0,NA)
astroCell_metadata$cellstatus <- ifelse(astroCell_metadata$AD_status>=0.25,1,astroCell_metadata$cellstatus)


AD_status_individual <- astroCell_metadata$AD_status_individual
cellstatus <- astroCell_metadata$cellstatus
plotData <- (table(AD_status_individual,cellstatus) / colSums(table(AD_status_individual,cellstatus)))|>as.data.frame()

medianscore <- lapply(unique(astroCell_metadata$projid),function(x){
    astroCell_metadata[astroCell_metadata$projid == x,"AD_status"] |>median(na.rm=T)
})|> unlist()

names(medianscore) <- unique(astroCell_metadata$projid)

astroCell_metadata$projid <- factor(astroCell_metadata$projid,levels = sort(medianscore) |>names())

## Supplementary figure 3
boxplotsAstro <- ggplot(astroCell_metadata,aes(as.factor(projid),AD_status, fill = AD_status_individual)) +
                        geom_boxplot(outlier.shape=NA) + geom_hline(yintercept = 0.25) +  geom_hline(yintercept = -0.25) + theme_minimal() + 
                        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
                        scale_fill_manual(values=c("blue", "red"))+ labs(y = "Cell projected AD score (astrocytes)", x = "377 donors sorted on median Cell projected AD score")

ggsave(boxplotsAstro, file = "Figures/Supplementary/Astrocyte_boxplots.pdf",width = 20, height = 7)


## Analysis
astro_cells <- raw_cells
astro_scores <- astroCell_metadata[astroCell_metadata$celltype %in% astro_cells@meta.data$cell_type_high_resolution,]
seurat_obj <- CreateSeuratObject(counts = astro_cells@assays$RNA@counts[,rownames(astro_scores)])
seurat_obj@meta.data$projid <- astro_scores$projid
seurat_obj@meta.data$AD_status_individual <- astro_scores$AD_status_individual
seurat_obj@meta.data$celltypes <- astro_scores$celltype
seurat_obj@meta.data$AD_status <- astro_scores$AD_status
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj@meta.data$cell_AD_status_binarized <- ifelse(seurat_obj@meta.data$AD_status<=-0.25,0,NA)
seurat_obj@meta.data$cell_AD_status_binarized <- ifelse(seurat_obj@meta.data$AD_status>=0.25,1,seurat_obj@meta.data$cell_AD_status_binarized)

## Filter only the healthy cells
seurat_obj_healthy_only <- subset(seurat_obj, subset = cell_AD_status_binarized == 0)

## Perfrom DEA between AD and non-AD, but only using their "healthy-like" cells
find <- FindMarkers(seurat_obj_healthy_only, group.by = "AD_status_individual",ident.1 = "0",logfc.threshold = 0.1)
table(find$p_val_adj < 0.01)

876 / 108

write.xlsx(find, file = "Results/DEA_Astrocytes_onDonorStatus_onlyHealthyCells(STable3).xlsx")

#Get cells that match the diagnosis of the individual
seurat_obj@meta.data$AD_status_individual <- as.numeric(seurat_obj@meta.data$AD_status_individual)
seurat_obj@meta.data$AD_status_individual <- seurat_obj@meta.data$AD_status_individual-1
seurat_obj@meta.data$matchCheck <- as.numeric(seurat_obj@meta.data$AD_status_individual) - seurat_obj@meta.data$cell_AD_status_binarized
## filter cells where diagnosis matches cell phenotype
seurat_obj_matching <- subset(seurat_obj, subset = matchCheck == 0)

#DGE between AD vs healthy where diagnosis matches cell phenotype
find_matching <- FindMarkers(seurat_obj_matching, group.by = "AD_status_individual",ident.1 = "0",logfc.threshold = 0.1)
table(find_matching$p_val_adj <= 0.01)

## supplementary table 2
write.xlsx(find_matching, file = "Results/DEA_Astrocytes_HealhydonorHealhyCells_vs_ADdonor_ADCells(STable4).xlsx")





