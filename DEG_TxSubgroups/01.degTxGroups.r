library(Seurat)
library(dplyr)
library(DESeq2)
library(ggpubr)
library(gridExtra)
library(ggrepel)


################################################################################
############################# FUNCTIONS ########################################
################################################################################

#' Safely run a DESeq2 analysis and return results
#'
#' This function wraps the DESeq2 workflow in a tryCatch block. If the analysis
#' is successful, it returns the results table. If it fails, it returns NULL
#' and prints an informative warning.
#'
#' @param count_data The raw count matrix.
#' @param col_data The column data (metadata) for the samples.
#' @param design_formula The design formula for the model.
#' @param contrast_name A character string name for the contrast, used for logging.
#' @return A DESeqResults object with added columns, or NULL if an error occurs.
run_deseq_safely <- function(count_data, col_data, design_formula, contrast_name) {
  # Use tryCatch to handle potential errors during the DESeq run
  res <- tryCatch({
    # 1. Create DESeqDataSet
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data,
                                          colData = col_data,
                                          design = design_formula)
    # 2. Run DESeq
    dds <- DESeq2::DESeq(dds)
    # 3. Get results
    results_table <- DESeq2::results(dds)
    # 4. Add custom columns
    results_table$minusLogPAdj <- -log10(results_table$padj)
    results_table$geneName <- rownames(results_table)
    # Print a success message
    message(paste("Successfully completed DESeq2 analysis for contrast:", contrast_name))
    return(results_table)
  }, error = function(e) {
    # This block executes if an error occurs in the 'try' block
    warning(paste("--> SKIPPED:", contrast_name, "contrast for this cell type due to an error.\n    Error message:", e$message))
    # Return NULL to indicate failure
    return(NULL)
  })
  return(res)
}

degPipeline <- function(cell_type, seurat_obj, merged_df_final, min_cells_per_donor = 30) {
  # Initialize an empty list to store successful results
  results_list <- list()
  #---- 1. Subgroups Contrast ----
  seurat_obj_currCellType_subgroups <- subset_seurat_subgroups(seurat_obj, cell_type, merged_df_final)
  pb_subgroups <- prepare_pseudobulk(seurat_obj_currCellType_subgroups, merged_df_final)
  pb_subgroups$coldata$subgroup <- factor(pb_subgroups$coldata$subgroup)
  pb_subgroups$coldata$subgroup <- relevel(pb_subgroups$coldata$subgroup , ref="PathNonAD-TxNonAD")
  res_subgroups <- run_deseq_safely(
    count_data = pb_subgroups$counts_donor,
    col_data = pb_subgroups$coldata,
    design_formula = ~ pmi + msex + age_death + batch + APOE_status + subgroup,
    contrast_name = "subgroups"
  )
  # Only add to the list if the result is not NULL
  if (!is.null(res_subgroups)) {
    results_list$subgroups <- res_subgroups
  }
  rm(seurat_obj_currCellType_subgroups, pb_subgroups)
  gc()
  #---- 2. Case/Control Contrast ----
  seurat_obj_currCellType_caseControl <- subset_seurat_caseControl(seurat_obj, cell_type, merged_df_final)
  pb_casecontrol <- prepare_pseudobulk(seurat_obj_currCellType_caseControl, merged_df_final)
  pb_casecontrol$coldata$AD_status <- factor(pb_casecontrol$coldata$AD_status)
  pb_casecontrol$coldata$AD_status <- relevel(pb_casecontrol$coldata$AD_status, ref = "0")
  res_casecontrol <- run_deseq_safely(
    count_data = pb_casecontrol$counts_donor,
    col_data = pb_casecontrol$coldata,
    design_formula = ~ pmi + msex + age_death + batch + APOE_status + AD_status,
    contrast_name = "casecontrol"
  )
  if (!is.null(res_casecontrol)) {
    results_list$casecontrol <- res_casecontrol
  }
  rm(seurat_obj_currCellType_caseControl, pb_casecontrol)
  gc()
  #---- 3. Cluster Contrast ----
  seurat_obj_currCellType_cluster <- subset_seurat_cluster(seurat_obj, cell_type, merged_df_final)
  pb_cluster <- prepare_pseudobulk(seurat_obj_currCellType_cluster, merged_df_final)
  pb_cluster$coldata$cluster_labels_for_big_heatmap <- factor(pb_cluster$coldata$cluster_labels_for_big_heatmap)
  pb_cluster$coldata$cluster_labels_for_big_heatmap <- relevel(pb_cluster$coldata$cluster_labels_for_big_heatmap, ref = "K1")
  res_cluster <- run_deseq_safely(
    count_data = pb_cluster$counts_donor,
    col_data = pb_cluster$coldata,
    design_formula = ~ pmi + msex + age_death + batch + APOE_status + cluster_labels_for_big_heatmap,
    contrast_name = "cluster"
  )
  if (!is.null(res_cluster)) {
    results_list$cluster <- res_cluster
  }
  rm(seurat_obj_currCellType_cluster, pb_cluster)
  gc()
  #---- 4. Shuffled Contrast ----
  seurat_obj_currCellType_shuf <- subset_seurat_shuf(seurat_obj, cell_type, merged_df_final)
  pb_shuf <- prepare_pseudobulk(seurat_obj_currCellType_shuf, merged_df_final)
  pb_shuf$coldata$shuf_group <- factor(pb_shuf$coldata$shuf_group)
  pb_shuf$coldata$shuf_group <- relevel(pb_shuf$coldata$shuf_group, ref = "GroupA")
  res_shuf <- run_deseq_safely(
    count_data = pb_shuf$counts_donor,
    col_data = pb_shuf$coldata,
    design_formula = ~ pmi + msex + age_death + batch + APOE_status + shuf_group,
    contrast_name = "shuf"
  )
  if (!is.null(res_shuf)) {
    results_list$shuf <- res_shuf
  }
  rm(seurat_obj_currCellType_shuf, pb_shuf)
  gc()
  #---- 5. Discordant Contrast ----
  seurat_obj_currCellType_discordant <- subset_seurat_discordant(seurat_obj, cell_type, merged_df_final)
  pb_discordant <- prepare_pseudobulk(seurat_obj_currCellType_discordant, merged_df_final)
  pb_discordant$coldata$subgroup <- factor(pb_discordant$coldata$subgroup)
  pb_discordant$coldata$subgroup <- relevel(pb_discordant$coldata$subgroup, ref = "PathNonAD-TxAD")
  res_discordant <- run_deseq_safely(
    count_data = pb_discordant$counts_donor,
    col_data = pb_discordant$coldata,
    design_formula = ~ pmi + msex + age_death + batch + APOE_status + subgroup,
    contrast_name = "discordant"
  )
  if (!is.null(res_discordant)) {
    results_list$discordant <- res_discordant
  }
  return(results_list)
}

degPipelineNoApoeCovariate <- function(cell_type, seurat_obj, merged_df_final, min_cells_per_donor = 30) {
  # Initialize an empty list to store successful results
  results_list <- list()
  #---- 1. Subgroups Contrast ----
  seurat_obj_currCellType_subgroups <- subset_seurat_subgroups(seurat_obj, cell_type, merged_df_final)
  pb_subgroups <- prepare_pseudobulk(seurat_obj_currCellType_subgroups, merged_df_final)
  pb_subgroups$coldata$subgroup <- factor(pb_subgroups$coldata$subgroup)
  pb_subgroups$coldata$subgroup <- relevel(pb_subgroups$coldata$subgroup , ref="PathNonAD-TxNonAD")
  res_subgroups <- run_deseq_safely(
    count_data = pb_subgroups$counts_donor,
    col_data = pb_subgroups$coldata,
    design_formula = ~ pmi + msex + age_death + batch + subgroup,
    contrast_name = "subgroups"
  )
  # Only add to the list if the result is not NULL
  if (!is.null(res_subgroups)) {
    results_list$subgroups <- res_subgroups
  }
  rm(seurat_obj_currCellType_subgroups, pb_subgroups)
  gc()
  #---- 2. Case/Control Contrast ----
  seurat_obj_currCellType_caseControl <- subset_seurat_caseControl(seurat_obj, cell_type, merged_df_final)
  pb_casecontrol <- prepare_pseudobulk(seurat_obj_currCellType_caseControl, merged_df_final)
  pb_casecontrol$coldata$AD_status <- factor(pb_casecontrol$coldata$AD_status)
  pb_casecontrol$coldata$AD_status <- relevel(pb_casecontrol$coldata$AD_status, ref = "0")
  res_casecontrol <- run_deseq_safely(
    count_data = pb_casecontrol$counts_donor,
    col_data = pb_casecontrol$coldata,
    design_formula = ~ pmi + msex + age_death + batch + AD_status,
    contrast_name = "casecontrol"
  )
  if (!is.null(res_casecontrol)) {
    results_list$casecontrol <- res_casecontrol
  }
  rm(seurat_obj_currCellType_caseControl, pb_casecontrol)
  gc()
  #---- 3. Cluster Contrast ----
  seurat_obj_currCellType_cluster <- subset_seurat_cluster(seurat_obj, cell_type, merged_df_final)
  pb_cluster <- prepare_pseudobulk(seurat_obj_currCellType_cluster, merged_df_final)
  pb_cluster$coldata$cluster_labels_for_big_heatmap <- factor(pb_cluster$coldata$cluster_labels_for_big_heatmap)
  pb_cluster$coldata$cluster_labels_for_big_heatmap <- relevel(pb_cluster$coldata$cluster_labels_for_big_heatmap, ref = "K1")
  res_cluster <- run_deseq_safely(
    count_data = pb_cluster$counts_donor,
    col_data = pb_cluster$coldata,
    design_formula = ~ pmi + msex + age_death + batch + cluster_labels_for_big_heatmap,
    contrast_name = "cluster"
  )
  if (!is.null(res_cluster)) {
    results_list$cluster <- res_cluster
  }
  rm(seurat_obj_currCellType_cluster, pb_cluster)
  gc()
  #---- 4. Shuffled Contrast ----
  seurat_obj_currCellType_shuf <- subset_seurat_shuf(seurat_obj, cell_type, merged_df_final)
  pb_shuf <- prepare_pseudobulk(seurat_obj_currCellType_shuf, merged_df_final)
  pb_shuf$coldata$shuf_group <- factor(pb_shuf$coldata$shuf_group)
  pb_shuf$coldata$shuf_group <- relevel(pb_shuf$coldata$shuf_group, ref = "GroupA")
  res_shuf <- run_deseq_safely(
    count_data = pb_shuf$counts_donor,
    col_data = pb_shuf$coldata,
    design_formula = ~ pmi + msex + age_death + batch + shuf_group,
    contrast_name = "shuf"
  )
  if (!is.null(res_shuf)) {
    results_list$shuf <- res_shuf
  }
  rm(seurat_obj_currCellType_shuf, pb_shuf)
  gc()
  #---- 5. Discordant Contrast ----
  seurat_obj_currCellType_discordant <- subset_seurat_discordant(seurat_obj, cell_type, merged_df_final)
  pb_discordant <- prepare_pseudobulk(seurat_obj_currCellType_discordant, merged_df_final)
  pb_discordant$coldata$subgroup <- factor(pb_discordant$coldata$subgroup)
  pb_discordant$coldata$subgroup <- relevel(pb_discordant$coldata$subgroup, ref = "PathNonAD-TxAD")
  res_discordant <- run_deseq_safely(
    count_data = pb_discordant$counts_donor,
    col_data = pb_discordant$coldata,
    design_formula = ~ pmi + msex + age_death + batch + subgroup,
    contrast_name = "discordant"
  )
  if (!is.null(res_discordant)) {
    results_list$discordant <- res_discordant
  }
  return(results_list)
}


degPipelineNoApoeNoBatchCovariate <- function(cell_type, seurat_obj, merged_df_final, min_cells_per_donor = 30) {
  # Initialize an empty list to store successful results
  results_list <- list()
  #---- 1. Subgroups Contrast ----
  seurat_obj_currCellType_subgroups <- subset_seurat_subgroups(seurat_obj, cell_type, merged_df_final)
  pb_subgroups <- prepare_pseudobulk(seurat_obj_currCellType_subgroups, merged_df_final)
  pb_subgroups$coldata$subgroup <- factor(pb_subgroups$coldata$subgroup)
  pb_subgroups$coldata$subgroup <- relevel(pb_subgroups$coldata$subgroup , ref="PathNonAD-TxNonAD")
  res_subgroups <- run_deseq_safely(
    count_data = pb_subgroups$counts_donor,
    col_data = pb_subgroups$coldata,
    design_formula = ~ pmi + msex + age_death + subgroup,
    contrast_name = "subgroups"
  )
  # Only add to the list if the result is not NULL
  if (!is.null(res_subgroups)) {
    results_list$subgroups <- res_subgroups
  }
  rm(seurat_obj_currCellType_subgroups, pb_subgroups)
  gc()
  #---- 2. Case/Control Contrast ----
  seurat_obj_currCellType_caseControl <- subset_seurat_caseControl(seurat_obj, cell_type, merged_df_final)
  pb_casecontrol <- prepare_pseudobulk(seurat_obj_currCellType_caseControl, merged_df_final)
  pb_casecontrol$coldata$AD_status <- factor(pb_casecontrol$coldata$AD_status)
  pb_casecontrol$coldata$AD_status <- relevel(pb_casecontrol$coldata$AD_status, ref = "0")
  res_casecontrol <- run_deseq_safely(
    count_data = pb_casecontrol$counts_donor,
    col_data = pb_casecontrol$coldata,
    design_formula = ~ pmi + msex + age_death + AD_status,
    contrast_name = "casecontrol"
  )
  if (!is.null(res_casecontrol)) {
    results_list$casecontrol <- res_casecontrol
  }
  rm(seurat_obj_currCellType_caseControl, pb_casecontrol)
  gc()
  #---- 3. Cluster Contrast ----
  seurat_obj_currCellType_cluster <- subset_seurat_cluster(seurat_obj, cell_type, merged_df_final)
  pb_cluster <- prepare_pseudobulk(seurat_obj_currCellType_cluster, merged_df_final)
  pb_cluster$coldata$cluster_labels_for_big_heatmap <- factor(pb_cluster$coldata$cluster_labels_for_big_heatmap)
  pb_cluster$coldata$cluster_labels_for_big_heatmap <- relevel(pb_cluster$coldata$cluster_labels_for_big_heatmap, ref = "K1")
  res_cluster <- run_deseq_safely(
    count_data = pb_cluster$counts_donor,
    col_data = pb_cluster$coldata,
    design_formula = ~ pmi + msex + age_death + cluster_labels_for_big_heatmap,
    contrast_name = "cluster"
  )
  if (!is.null(res_cluster)) {
    results_list$cluster <- res_cluster
  }
  rm(seurat_obj_currCellType_cluster, pb_cluster)
  gc()
  #---- 4. Shuffled Contrast ----
  seurat_obj_currCellType_shuf <- subset_seurat_shuf(seurat_obj, cell_type, merged_df_final)
  pb_shuf <- prepare_pseudobulk(seurat_obj_currCellType_shuf, merged_df_final)
  pb_shuf$coldata$shuf_group <- factor(pb_shuf$coldata$shuf_group)
  pb_shuf$coldata$shuf_group <- relevel(pb_shuf$coldata$shuf_group, ref = "GroupA")
  res_shuf <- run_deseq_safely(
    count_data = pb_shuf$counts_donor,
    col_data = pb_shuf$coldata,
    design_formula = ~ pmi + msex + age_death + shuf_group,
    contrast_name = "shuf"
  )
  if (!is.null(res_shuf)) {
    results_list$shuf <- res_shuf
  }
  rm(seurat_obj_currCellType_shuf, pb_shuf)
  gc()
  #---- 5. Discordant Contrast ----
  seurat_obj_currCellType_discordant <- subset_seurat_discordant(seurat_obj, cell_type, merged_df_final)
  pb_discordant <- prepare_pseudobulk(seurat_obj_currCellType_discordant, merged_df_final)
  pb_discordant$coldata$subgroup <- factor(pb_discordant$coldata$subgroup)
  pb_discordant$coldata$subgroup <- relevel(pb_discordant$coldata$subgroup, ref = "PathNonAD-TxAD")
  res_discordant <- run_deseq_safely(
    count_data = pb_discordant$counts_donor,
    col_data = pb_discordant$coldata,
    design_formula = ~ pmi + msex + age_death + subgroup,
    contrast_name = "discordant"
  )
  if (!is.null(res_discordant)) {
    results_list$discordant <- res_discordant
  }
  return(results_list)
}

subset_seurat_subgroups <- function(seurat_obj, cell_type, merged_df_final) {
  currSeurat <- subset(seurat_obj, cell_type_high_resolution == cell_type)
  currSeurat$cell <- rownames(currSeurat@meta.data)
  currSeurat@meta.data$projid <- as.integer(currSeurat@meta.data$projid)
  currSeurat <- subset(
    currSeurat,
    cells = rownames(currSeurat@meta.data)[
      currSeurat@meta.data$projid %in% merged_df_final$projid
    ]
  )
  currSeurat@meta.data <- currSeurat@meta.data %>%
    left_join(merged_df_final, by = "projid")
  rownames(currSeurat@meta.data) <- currSeurat@meta.data$cell
  currSeurat <- subset(currSeurat, subset = subgroup %in% c("PathAD-TxAD", "PathNonAD-TxNonAD"))
  return(currSeurat)
}

subset_seurat_discordant <- function(seurat_obj, cell_type, merged_df_final) {
  currSeurat <- subset(seurat_obj, cell_type_high_resolution == cell_type)
  currSeurat$cell <- rownames(currSeurat@meta.data)
  currSeurat@meta.data$projid <- as.integer(currSeurat@meta.data$projid)
  currSeurat <- subset(
    currSeurat,
    cells = rownames(currSeurat@meta.data)[
      currSeurat@meta.data$projid %in% merged_df_final$projid
    ]
  )
  currSeurat@meta.data <- currSeurat@meta.data %>%
    left_join(merged_df_final, by = "projid")
  rownames(currSeurat@meta.data) <- currSeurat@meta.data$cell
  currSeurat <- subset(currSeurat, subset = subgroup %in% c("PathAD-TxNonAD", "PathNonAD-TxAD"))
  return(currSeurat)
}

subset_seurat_caseControl <- function(seurat_obj, cell_type, merged_df_final) {
  currSeurat <- subset(seurat_obj, cell_type_high_resolution == cell_type)
  currSeurat$cell <- rownames(currSeurat@meta.data)
  currSeurat@meta.data$projid <- as.integer(currSeurat@meta.data$projid)
  currSeurat <- subset(
    currSeurat,
    cells = rownames(currSeurat@meta.data)[
      currSeurat@meta.data$projid %in% merged_df_final$projid
    ]
  )
  currSeurat@meta.data <- currSeurat@meta.data %>%
    left_join(merged_df_final, by = "projid")
  rownames(currSeurat@meta.data) <- currSeurat@meta.data$cell
  #all donors are either case or control, no further subsetting required
  return(currSeurat)
}

subset_seurat_cluster <- function(seurat_obj, cell_type, merged_df_final) {
  currSeurat <- subset(seurat_obj, cell_type_high_resolution == cell_type)
  currSeurat$cell <- rownames(currSeurat@meta.data)
  currSeurat@meta.data$projid <- as.integer(currSeurat@meta.data$projid)
  currSeurat <- subset(
    currSeurat,
    cells = rownames(currSeurat@meta.data)[
      currSeurat@meta.data$projid %in% merged_df_final$projid
    ]
  )
  currSeurat@meta.data <- currSeurat@meta.data %>%
    left_join(merged_df_final, by = "projid")
  rownames(currSeurat@meta.data) <- currSeurat@meta.data$cell
  currSeurat <- subset(currSeurat, subset = cluster_labels_for_big_heatmap %in% c("K1", "K3"))
  return(currSeurat)
}

subset_seurat_shuf <- function(seurat_obj, cell_type, maerged_df_final) {
  currSeurat <- subset(seurat_obj, cell_type_high_resolution == cell_type)
  currSeurat$cell <- rownames(currSeurat@meta.data)
  currSeurat@meta.data$projid <- as.integer(currSeurat@meta.data$projid)
  currSeurat <- subset(
    currSeurat,
    cells = rownames(currSeurat@meta.data)[
      currSeurat@meta.data$projid %in% merged_df_final$projid
    ]
  )
  currSeurat@meta.data <- currSeurat@meta.data %>%
    left_join(merged_df_final, by = "projid")
  rownames(currSeurat@meta.data) <- currSeurat@meta.data$cell
  #all donors are either GroupA or GroupB, no further subsetting required
  return(currSeurat)
}

prepare_pseudobulk <- function(obj, merged_df_final, min_cells_per_donor = 30) {
  cell_counts_per_donor <- table(obj@meta.data$projid)
  counts <- obj@assays$RNA@counts
  cell_to_donor <- obj@meta.data$projid
  names(cell_to_donor) <- colnames(obj)
  donors <- unique(cell_to_donor)
  counts_donor <- matrix(0, nrow = nrow(counts), ncol = length(donors),
                         dimnames = list(rownames(counts), donors))
  for (i in seq_along(donors)) {
    donor_cells <- which(cell_to_donor == donors[i])
    sub_counts <- counts[, donor_cells, drop = FALSE]
    counts_donor[, i] <- rowSums(sub_counts)
  }
  donors_keep <- names(cell_counts_per_donor)[cell_counts_per_donor >= min_cells_per_donor]
  counts_donor <- counts_donor[, donors_keep, drop = FALSE]
  coldata <- merged_df_final %>% filter(projid %in% donors_keep)
  coldata <- coldata[match(colnames(counts_donor), coldata$projid), ]
  return(list(counts_donor = counts_donor, coldata = coldata))
}

################################################################################
########################## START ANALYSIS ######################################
################################################################################

set.seed(1612)
merged_df_final <- readRDS("merged_df_final_11_03.rds")
merged_df_final$projid <- as.integer(merged_df_final$projid)
merged_df_final$shuf_group <- sample(c("GroupA", "GroupB"),
                                     size = nrow(merged_df_final),
                                     replace = TRUE)

proj_batch <- read.csv("meta_data_with_batch.csv")
merged_df_final <- merged_df_final %>%
  left_join(proj_batch %>% select(projid, batch), 
            by = "projid")

saveRDS(merged_df_final, "merged_df_final_12_18.rds")
merged_df_final <- readRDS("merged_df_final_12_18.rds")

immuneSeurat <- readRDS("../Data/mathys/Immune_cells.rds") #raw seurat object
micResults <- degPipeline("Mic P2RY12", immuneSeurat, merged_df_final)
saveRDS(micResults, "Results/Mic P2RY12.rds")
rm(immuneSeurat)
gc()

exc1Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set1.rds")
ExcL23CBLN2LINC02306_Results <- degPipeline("Exc L2-3 CBLN2 LINC02306", exc1Seurat, merged_df_final)
saveRDS(ExcL23CBLN2LINC02306_Results, "Results/Exc L2-3 CBLN2 LINC02306.rds")
rm(exc1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
rorbSeurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L3-4 RORB CUX2")
rm(exc2Seurat)
gc()
Exc.L3.4.RORB.CUX2_deg_results <- degPipeline("Exc L3-4 RORB CUX2", rorbSeurat, merged_df_final)
saveRDS(Exc.L3.4.RORB.CUX2_deg_results, "Results/Exc L3-4 RORB CUX2.rds")
rm(rorbSeurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
GABRG1Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L4-5 RORB GABRG1")
rm(exc2Seurat)
gc()
Exc.L4.5.RORB.GABRG1_deg_results <- degPipeline("Exc L4-5 RORB GABRG1", GABRG1Seurat, merged_df_final)
saveRDS(Exc.L4.5.RORB.GABRG1_deg_results, "Results/Exc L4-5 RORB GABRG1.rds")
rm(GABRG1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
PLCH1Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L3-5 RORB PLCH1")
rm(exc2Seurat)
gc()
Exc.L3.5.RORB.PLCH1_deg_results <- degPipeline("Exc L3-5 RORB PLCH1", PLCH1Seurat, merged_df_final)
saveRDS(Exc.L3.5.RORB.PLCH1_deg_results, "Results/Exc L3-5 RORB PLCH1.rds")
rm(PLCH1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
IL1RAPL2Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L4-5 RORB IL1RAPL2")
rm(exc2Seurat)
gc()
Exc.L4.5.RORB.IL1RAPL2_deg_results <- degPipeline("Exc L4-5 RORB IL1RAPL2", IL1RAPL2Seurat, merged_df_final)
saveRDS(Exc.L4.5.RORB.IL1RAPL2_deg_results, "Results/Exc L4-5 RORB IL1RAPL2.rds")
rm(IL1RAPL2Seurat)
gc()

exc3Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set3.rds")

## too few cells, matrix not full rank
#Exc.L5.ET_deg_results <- degPipeline("Exc L5 ET", exc3Seurat, merged_df_final)
#saveRDS(Exc.L5.ET_deg_results, "Results/Exc L5 ET.rds")
Exc.L5.6.RORB.LINC02196_deg_results <- degPipeline("Exc L5-6 RORB LINC02196", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.RORB.LINC02196_deg_results, "Results/Exc L5-6 RORB LINC02196.rds")
Exc.L5.6.IT.Car3_deg_results <- degPipeline("Exc L5/6 IT Car3", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.IT.Car3_deg_results, "Results/Exc L5-6 IT Car3.rds")
Exc.L5.6.NP_deg_results <- degPipeline("Exc L5/6 NP", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.NP_deg_results, "Results/Exc L5-6 NP.rds")
Exc.L6.CT_deg_results <- degPipeline("Exc L6 CT", exc3Seurat, merged_df_final)
saveRDS(Exc.L6.CT_deg_results, "Results/Exc L6 CT.rds")
Exc.L6.THEMIS.NFIA_deg_results <- degPipeline("Exc L6 THEMIS NFIA", exc3Seurat, merged_df_final)
saveRDS(Exc.L6.THEMIS.NFIA_deg_results, "Results/Exc L6 THEMIS NFIA.rds")
Exc.L6b_deg_results <- degPipeline("Exc L6b", exc3Seurat, merged_df_final)
saveRDS(Exc.L6b_deg_results, "Results/Exc L6b.rds")
Exc.NRGN_deg_results <- degPipeline("Exc NRGN", exc3Seurat, merged_df_final)
saveRDS(Exc.NRGN_deg_results, "Results/Exc NRGN.rds")
Exc.RELN.CHD7_deg_results <- degPipeline("Exc RELN CHD7", exc3Seurat, merged_df_final)
saveRDS(Exc.RELN.CHD7_deg_results, "Results/Exc RELN CHD7.rds")
rm(exc3Seurat)
gc()

inh_list <- c(
  "Inh PVALB CA8 (Chandelier)",
  "Inh CUX2 MSR1",
  "Inh RYR3 TSHZ2",
  "Inh L3-5 SST MAFB",
  "Inh PVALB HTR4",
  "Inh ENOX2 SPHKAP",
  ##"Inh LAMP5 RELN", # too few cells to run
  "Inh VIP CLSTN2",
  ##"Inh GPC5 RIT2", # too few cells to run
  "Inh VIP ABI3BP",
  ##"Inh L1 PAX6 CA4", # too few cells to run
  "Inh PVALB SULF1",
  "Inh ALCAM TRPM3",
  ##"Inh L5-6 SST TH", # too few cells to run
  "Inh PTPRK FAM19A1",
  "Inh LAMP5 NRG1 (Rosehip)",
  "Inh L1-6 LAMP5 CA13",
  "Inh VIP TSHZ2",
  "Inh SORCS1 TTN",
  ##"Inh FBN2 EPB41L4A", # too few cells to run
  ##"Inh L1-2 PAX6 SCGN",# too few cells to run
  ##"Inh L5-6 PVALB STON2",# too few cells to run
  "Inh VIP THSD7B"
  ##"Inh SGCD PDE3A" # too few cells to run
)

inhSeurat <- readRDS("../Data/mathys/Inhibitory_neurons.rds")
inh_results <- lapply(inh_list, function(m) {
  res <- degPipeline(m, inhSeurat, merged_df_final)
  file_name <- paste0("Results/", m, ".rds")
  saveRDS(res, file_name)
  return(res)
})
rm(inhSeurat)
gc()

astSeurat <- readRDS("../Data/mathys/Astrocytes.rds")
astSeuratV5 <- UpdateSeuratObject(astSeurat)
Ast.GRM3_deg_results <- degPipeline("Ast GRM3", astSeuratV5, merged_df_final)
saveRDS(Ast.GRM3_deg_results, "Results/Ast GRM3.rds")
Ast.CHI3L1_deg_results <- degPipeline("Ast CHI3L1", astSeuratV5, merged_df_final)
saveRDS(Ast.CHI3L1_deg_results, "Results/Ast CHI3L1.rds")
Ast.DPP10_deg_results <- degPipeline("Ast DPP10", astSeuratV5, merged_df_final)
saveRDS(Ast.DPP10_deg_results, "Results/Ast DPP10.rds")
rm(astSeurat)
rm(astSeuratV5)
gc()

oliSeurat <- readRDS("../Data/mathys/Oligodendrocytes.rds")
oliSeuratV5 <- UpdateSeuratObject(oliSeurat)
Oli_deg_results <- degPipeline("Oli", oliSeuratV5, merged_df_final)
saveRDS(Oli_deg_results, "Results/Oli.rds")
rm(oliSeurat)
rm(oliSeuratV5)
gc()

opcSeurat <- readRDS("../Data/mathys/OPCs.rds")
opc_deg_results <- degPipeline("OPC", opcSeurat, merged_df_final)
saveRDS(opc_deg_results, "Results/OPC.rds")
rm(opcSeurat)
gc()

vascSeurat <- readRDS("../Data/mathys/Vasculature_cells.rds")
end_deg_results <- degPipeline("End", vascSeurat, merged_df_final)# too few to run
saveRDS(end_deg_results, "Results/End.rds") # too few to run
Fib.FLRT2_deg_results <- degPipeline("Fib FLRT2", vascSeurat, merged_df_final)
saveRDS(Fib.FLRT2_deg_results, "Results/Fib FLRT2.rds")


################################################################################
######################## ANALYSIS w/o APOE #####################################
################################################################################

merged_df_final <- readRDS("merged_df_final_12_18.rds")

immuneSeurat <- readRDS("../Data/mathys/Immune_cells.rds") #raw seurat object
micResults <- degPipelineNoApoeCovariate("Mic P2RY12", immuneSeurat, merged_df_final)
saveRDS(micResults, "ResultsNoApoeCovariate/Mic P2RY12.rds")
rm(immuneSeurat)
gc()

exc1Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set1.rds")
ExcL23CBLN2LINC02306_Results <- degPipelineNoApoeCovariate("Exc L2-3 CBLN2 LINC02306", exc1Seurat, merged_df_final)
saveRDS(ExcL23CBLN2LINC02306_Results, "ResultsNoApoeCovariate/Exc L2-3 CBLN2 LINC02306.rds")
rm(exc1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
rorbSeurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L3-4 RORB CUX2")
rm(exc2Seurat)
gc()
Exc.L3.4.RORB.CUX2_deg_results <- degPipelineNoApoeCovariate("Exc L3-4 RORB CUX2", rorbSeurat, merged_df_final)
saveRDS(Exc.L3.4.RORB.CUX2_deg_results, "ResultsNoApoeCovariate/Exc L3-4 RORB CUX2.rds")
rm(rorbSeurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
GABRG1Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L4-5 RORB GABRG1")
rm(exc2Seurat)
gc()
Exc.L4.5.RORB.GABRG1_deg_results <- degPipelineNoApoeCovariate("Exc L4-5 RORB GABRG1", GABRG1Seurat, merged_df_final)
saveRDS(Exc.L4.5.RORB.GABRG1_deg_results, "ResultsNoApoeCovariate/Exc L4-5 RORB GABRG1.rds")
rm(GABRG1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
PLCH1Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L3-5 RORB PLCH1")
rm(exc2Seurat)
gc()
Exc.L3.5.RORB.PLCH1_deg_results <- degPipelineNoApoeCovariate("Exc L3-5 RORB PLCH1", PLCH1Seurat, merged_df_final)
saveRDS(Exc.L3.5.RORB.PLCH1_deg_results, "ResultsNoApoeCovariate/Exc L3-5 RORB PLCH1.rds")
rm(PLCH1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
IL1RAPL2Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L4-5 RORB IL1RAPL2")
rm(exc2Seurat)
gc()
Exc.L4.5.RORB.IL1RAPL2_deg_results <- degPipelineNoApoeCovariate("Exc L4-5 RORB IL1RAPL2", IL1RAPL2Seurat, merged_df_final)
saveRDS(Exc.L4.5.RORB.IL1RAPL2_deg_results, "ResultsNoApoeCovariate/Exc L4-5 RORB IL1RAPL2.rds")
rm(IL1RAPL2Seurat)
gc()

exc3Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set3.rds")

## too few cells, matrix not full rank
#Exc.L5.ET_deg_results <- degPipelineNoApoeCovariate("Exc L5 ET", exc3Seurat, merged_df_final)
#saveRDS(Exc.L5.ET_deg_results, "ResultsNoApoeCovariate/Exc L5 ET.rds")
Exc.L5.6.RORB.LINC02196_deg_results <- degPipelineNoApoeCovariate("Exc L5-6 RORB LINC02196", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.RORB.LINC02196_deg_results, "ResultsNoApoeCovariate/Exc L5-6 RORB LINC02196.rds")
Exc.L5.6.IT.Car3_deg_results <- degPipelineNoApoeCovariate("Exc L5/6 IT Car3", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.IT.Car3_deg_results, "ResultsNoApoeCovariate/Exc L5-6 IT Car3.rds")
Exc.L5.6.NP_deg_results <- degPipelineNoApoeCovariate("Exc L5/6 NP", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.NP_deg_results, "ResultsNoApoeCovariate/Exc L5-6 NP.rds")
Exc.L6.CT_deg_results <- degPipelineNoApoeCovariate("Exc L6 CT", exc3Seurat, merged_df_final)
saveRDS(Exc.L6.CT_deg_results, "ResultsNoApoeCovariate/Exc L6 CT.rds")
Exc.L6.THEMIS.NFIA_deg_results <- degPipelineNoApoeCovariate("Exc L6 THEMIS NFIA", exc3Seurat, merged_df_final)
saveRDS(Exc.L6.THEMIS.NFIA_deg_results, "ResultsNoApoeCovariate/Exc L6 THEMIS NFIA.rds")
Exc.L6b_deg_results <- degPipelineNoApoeCovariate("Exc L6b", exc3Seurat, merged_df_final)
saveRDS(Exc.L6b_deg_results, "ResultsNoApoeCovariate/Exc L6b.rds")
Exc.NRGN_deg_results <- degPipelineNoApoeCovariate("Exc NRGN", exc3Seurat, merged_df_final)
saveRDS(Exc.NRGN_deg_results, "ResultsNoApoeCovariate/Exc NRGN.rds")
Exc.RELN.CHD7_deg_results <- degPipelineNoApoeCovariate("Exc RELN CHD7", exc3Seurat, merged_df_final)
saveRDS(Exc.RELN.CHD7_deg_results, "ResultsNoApoeCovariate/Exc RELN CHD7.rds")
rm(exc3Seurat)
gc()


inhSeurat <- readRDS("../Data/mathys/Inhibitory_neurons.rds")
inh_results <- lapply(inh_list, function(m) {
  res <- degPipelineNoApoeCovariate(m, inhSeurat, merged_df_final)
  file_name <- paste0("ResultsNoApoeCovariate/", m, ".rds")
  saveRDS(res, file_name)
  return(res)
})
rm(inhSeurat)
gc()

astSeurat <- readRDS("../Data/mathys/Astrocytes.rds")
astSeuratV5 <- UpdateSeuratObject(astSeurat)
Ast.GRM3_deg_results <- degPipelineNoApoeCovariate("Ast GRM3", astSeuratV5, merged_df_final)
saveRDS(Ast.GRM3_deg_results, "ResultsNoApoeCovariate/Ast GRM3.rds")
Ast.CHI3L1_deg_results <- degPipelineNoApoeCovariate("Ast CHI3L1", astSeuratV5, merged_df_final)
saveRDS(Ast.CHI3L1_deg_results, "ResultsNoApoeCovariate/Ast CHI3L1.rds")
Ast.DPP10_deg_results <- degPipelineNoApoeCovariate("Ast DPP10", astSeuratV5, merged_df_final)
saveRDS(Ast.DPP10_deg_results, "ResultsNoApoeCovariate/Ast DPP10.rds")
rm(astSeurat)
rm(astSeuratV5)
gc()

oliSeurat <- readRDS("../Data/mathys/Oligodendrocytes.rds")
oliSeuratV5 <- UpdateSeuratObject(oliSeurat)
Oli_deg_results <- degPipelineNoApoeCovariate("Oli", oliSeuratV5, merged_df_final)
saveRDS(Oli_deg_results, "ResultsNoApoeCovariate/Oli.rds")
rm(oliSeurat)
rm(oliSeuratV5)
gc()

opcSeurat <- readRDS("../Data/mathys/OPCs.rds")
opc_deg_results <- degPipelineNoApoeCovariate("OPC", opcSeurat, merged_df_final)
saveRDS(opc_deg_results, "ResultsNoApoeCovariate/OPC.rds")
rm(opcSeurat)
gc()

vascSeurat <- readRDS("../Data/mathys/Vasculature_cells.rds")
#end_deg_results <- degPipelineNoApoeCovariate("End", vascSeurat, merged_df_final)# too few to run
#saveRDS(end_deg_results, "ResultsNoApoeCovariate/End.rds") # too few to run
Fib.FLRT2_deg_results <- degPipelineNoApoeCovariate("Fib FLRT2", vascSeurat, merged_df_final)
saveRDS(Fib.FLRT2_deg_results, "ResultsNoApoeCovariate/Fib FLRT2.rds")

################################################################################
################ ANALYSIS w/o APOE or Batch covariates #########################
################################################################################

merged_df_final <- readRDS("merged_df_final_12_18.rds")

immuneSeurat <- readRDS("../Data/mathys/Immune_cells.rds") #raw seurat object
micResults <- degPipelineNoApoeNoBatchCovariate("Mic P2RY12", immuneSeurat, merged_df_final)
saveRDS(micResults, "ResultsNoSeqBatchNoApoeCovariate/Mic P2RY12.rds")
rm(immuneSeurat)
gc()

exc1Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set1.rds")
ExcL23CBLN2LINC02306_Results <- degPipelineNoApoeNoBatchCovariate("Exc L2-3 CBLN2 LINC02306", exc1Seurat, merged_df_final)
saveRDS(ExcL23CBLN2LINC02306_Results, "ResultsNoSeqBatchNoApoeCovariate/Exc L2-3 CBLN2 LINC02306.rds")
rm(exc1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
rorbSeurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L3-4 RORB CUX2")
rm(exc2Seurat)
gc()
Exc.L3.4.RORB.CUX2_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L3-4 RORB CUX2", rorbSeurat, merged_df_final)
saveRDS(Exc.L3.4.RORB.CUX2_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L3-4 RORB CUX2.rds")
rm(rorbSeurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
GABRG1Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L4-5 RORB GABRG1")
rm(exc2Seurat)
gc()
Exc.L4.5.RORB.GABRG1_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L4-5 RORB GABRG1", GABRG1Seurat, merged_df_final)
saveRDS(Exc.L4.5.RORB.GABRG1_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L4-5 RORB GABRG1.rds")
rm(GABRG1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
PLCH1Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L3-5 RORB PLCH1")
rm(exc2Seurat)
gc()
Exc.L3.5.RORB.PLCH1_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L3-5 RORB PLCH1", PLCH1Seurat, merged_df_final)
saveRDS(Exc.L3.5.RORB.PLCH1_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L3-5 RORB PLCH1.rds")
rm(PLCH1Seurat)
gc()

exc2Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set2.rds")
IL1RAPL2Seurat <- subset(exc2Seurat, subset = cell_type_high_resolution == "Exc L4-5 RORB IL1RAPL2")
rm(exc2Seurat)
gc()
Exc.L4.5.RORB.IL1RAPL2_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L4-5 RORB IL1RAPL2", IL1RAPL2Seurat, merged_df_final)
saveRDS(Exc.L4.5.RORB.IL1RAPL2_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L4-5 RORB IL1RAPL2.rds")
rm(IL1RAPL2Seurat)
gc()

exc3Seurat <- readRDS("../Data/mathys/Excitatory_neurons_set3.rds")

## too few cells, matrix not full rank
#Exc.L5.ET_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L5 ET", exc3Seurat, merged_df_final)
#saveRDS(Exc.L5.ET_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L5 ET.rds")
Exc.L5.6.RORB.LINC02196_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L5-6 RORB LINC02196", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.RORB.LINC02196_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L5-6 RORB LINC02196.rds")
Exc.L5.6.IT.Car3_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L5/6 IT Car3", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.IT.Car3_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L5-6 IT Car3.rds")
Exc.L5.6.NP_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L5/6 NP", exc3Seurat, merged_df_final)
saveRDS(Exc.L5.6.NP_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L5-6 NP.rds")
Exc.L6.CT_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L6 CT", exc3Seurat, merged_df_final)
saveRDS(Exc.L6.CT_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L6 CT.rds")
Exc.L6.THEMIS.NFIA_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L6 THEMIS NFIA", exc3Seurat, merged_df_final)
saveRDS(Exc.L6.THEMIS.NFIA_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L6 THEMIS NFIA.rds")
Exc.L6b_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc L6b", exc3Seurat, merged_df_final)
saveRDS(Exc.L6b_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc L6b.rds")
Exc.NRGN_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc NRGN", exc3Seurat, merged_df_final)
saveRDS(Exc.NRGN_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc NRGN.rds")
Exc.RELN.CHD7_deg_results <- degPipelineNoApoeNoBatchCovariate("Exc RELN CHD7", exc3Seurat, merged_df_final)
saveRDS(Exc.RELN.CHD7_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Exc RELN CHD7.rds")
rm(exc3Seurat)
gc()

inh_list <- c(
  "Inh PVALB CA8 (Chandelier)",
  "Inh CUX2 MSR1",
  "Inh RYR3 TSHZ2",
  "Inh L3-5 SST MAFB",
  "Inh PVALB HTR4",
  "Inh ENOX2 SPHKAP",
  ##"Inh LAMP5 RELN", # too few cells to run
  "Inh VIP CLSTN2",
  ##"Inh GPC5 RIT2", # too few cells to run
  "Inh VIP ABI3BP",
  ##"Inh L1 PAX6 CA4", # too few cells to run
  "Inh PVALB SULF1",
  "Inh ALCAM TRPM3",
  ##"Inh L5-6 SST TH", # too few cells to run
  "Inh PTPRK FAM19A1",
  "Inh LAMP5 NRG1 (Rosehip)",
  "Inh L1-6 LAMP5 CA13",
  "Inh VIP TSHZ2",
  "Inh SORCS1 TTN",
  ##"Inh FBN2 EPB41L4A", # too few cells to run
  ##"Inh L1-2 PAX6 SCGN",# too few cells to run
  ##"Inh L5-6 PVALB STON2",# too few cells to run
  "Inh VIP THSD7B"
  ##"Inh SGCD PDE3A" # too few cells to run
)

inhSeurat <- readRDS("../Data/mathys/Inhibitory_neurons.rds")
inh_results <- lapply(inh_list, function(m) {
  res <- degPipelineNoApoeNoBatchCovariate(m, inhSeurat, merged_df_final)
  file_name <- paste0("ResultsNoSeqBatchNoApoeCovariate/", m, ".rds")
  saveRDS(res, file_name)
  return(res)
})
rm(inhSeurat)
gc()

astSeurat <- readRDS("../Data/mathys/Astrocytes.rds")
astSeuratV5 <- UpdateSeuratObject(astSeurat)
Ast.GRM3_deg_results <- degPipelineNoApoeNoBatchCovariate("Ast GRM3", astSeuratV5, merged_df_final)
saveRDS(Ast.GRM3_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Ast GRM3.rds")
Ast.CHI3L1_deg_results <- degPipelineNoApoeNoBatchCovariate("Ast CHI3L1", astSeuratV5, merged_df_final)
saveRDS(Ast.CHI3L1_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Ast CHI3L1.rds")
Ast.DPP10_deg_results <- degPipelineNoApoeNoBatchCovariate("Ast DPP10", astSeuratV5, merged_df_final)
saveRDS(Ast.DPP10_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Ast DPP10.rds")
rm(astSeurat)
rm(astSeuratV5)
gc()

oliSeurat <- readRDS("../Data/mathys/Oligodendrocytes.rds")
oliSeuratV5 <- UpdateSeuratObject(oliSeurat)
Oli_deg_results <- degPipelineNoApoeNoBatchCovariate("Oli", oliSeuratV5, merged_df_final)
saveRDS(Oli_deg_results, "ResultsNoSeqBatchNoApoeCovariate/Oli.rds")
rm(oliSeurat)
rm(oliSeuratV5)
gc()

opcSeurat <- readRDS("../Data/mathys/OPCs.rds")
opc_deg_results <- degPipelineNoApoeNoBatchCovariate("OPC", opcSeurat, merged_df_final)
saveRDS(opc_deg_results, "ResultsNoSeqBatchNoApoeCovariate/OPC.rds")
rm(opcSeurat)
gc()
