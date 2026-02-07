library(uwot)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)

################################################################################
############################### DF SETUP #######################################
################################################################################

merged_df_final <- readRDS("Data/merged_df_final_11_03.rds")
rownames(merged_df_final) <- merged_df_final$projid
cols_to_keep <- c(
  "projid", "cogn_global_lv", "age_death", "msex", "gpath", "plaq_n",
  "amyloid", "tangles", "AD_status", "CR", "APOE_status", "cogn_decline",  
  "cluster_labels_for_big_heatmap", "subgroup"
)
merged_df_clean <- merged_df_final[, cols_to_keep, drop = FALSE]
cell_scores_list <- load_celltype_revised_phenotype(c("AD_status", "amyloid", "tangles", "APOE_status", "plaq_n"), dataset = "mathys")
combine_cell_scores <- function(cell_scores_list) {
  # For each endophenotype: rename columns and return matrix
  renamed_mats <- lapply(names(cell_scores_list), function(pheno) {
    mat <- cell_scores_list[[pheno]]
    colnames(mat) <- paste("CPP", pheno, colnames(mat), sep = ".")
    mat
  })
  combined_mat <- do.call(cbind, renamed_mats)
  return(combined_mat)
}

combined_mat <- combine_cell_scores(cell_scores_list)
pca <- prcomp(combined_mat, scale. = TRUE)
num_pcs <- 20
df_pca <- data.frame(
  projid = rownames(combined_mat),
  pca$x[, 1:num_pcs],
  check.names = FALSE
)

colnames(df_pca)[-1] <- paste0("PC", seq_len(num_pcs))

combined_df <- data.frame(
  projid = rownames(combined_mat),
  combined_mat,
  check.names = FALSE
)

tmp <- merge(
  merged_df_clean,
  combined_df,
  by = "projid",
  all.x = TRUE,
  sort = FALSE
)

merged_annotated <- merge(
  tmp,
  df_pca,
  by = "projid",
  all.x = TRUE,
  sort = FALSE
)

set.seed(123)  # reproducibility
umap_res <- umap(combined_mat, n_neighbors = 30, min_dist = 0.3, metric = "euclidean")

df_umap <- data.frame(
  projid = rownames(combined_mat),
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2]
)

# --- 2. Merge with merged_annotated ---
merged_annotated <- merge(
  merged_annotated,
  df_umap,
  by = "projid",
  all.x = TRUE,
  sort = FALSE
)

################################################################################
############################### META VIZ #######################################
################################################################################
a <- ggscatter(merged_annotated, "PC1" ,"PC2", color="AD_status")+ scale_color_gradient(low="blue", high="red") + 
  ggtitle("Donor AD Status")
subgroup_cols <- c(
  "PathNonAD-TxMid" = "#fcbf49",
  "PathAD-TxAD" = "#ad2524",
  "PathAD-TxMid" = "#f38b20",
  "PathAD-TxNonAD" ="#90a955",
  "PathNonAD-TxNonAD" = "#264653",
  "PathNonAD-TxAD" = "#9d4edd"
)
b <- ggscatter(merged_annotated, "PC1", "PC2", color="subgroup", pal=subgroup_cols) + ggtitle("CPP-based Tx Subgroup")
cluster_cols <- c(
  "K1" = "#264653",
  "K2" = "#f38b20",
  "K3" = "#ad2524"
)
c <- ggscatter(merged_annotated, "PC1", "PC2", color="cluster_labels_for_big_heatmap", pal=cluster_cols) + ggtitle("CPP-based Tx Cluster")
d <- ggscatter(merged_annotated, "PC1", "PC2", color="tangles") + scale_color_gradient(low="gray", high="orange")
e <- ggscatter(merged_annotated, "PC1", "PC2", color="gpath") + scale_color_gradient(low="gray", high="darkgreen")
f <- ggscatter(merged_annotated, "PC1", "PC2", color="cogn_decline") + scale_color_gradient(low="gray", high = "purple")
g <- ggscatter(merged_annotated, "PC1", "PC2", color="age_death") + scale_color_gradient(low="gray", high = "brown")
h <- ggscatter(merged_annotated, "PC1", "PC2", color="msex") + scale_color_gradient(low="#4895ef", high = "#fb6f92")

pdf("pdf/subtyping/pcaWithLegends.pdf", height=12, width=15)
grid.arrange(a,b,c,d,e,f,g,h,nrow=3)
dev.off()

a <- a + theme(legend.position="none") + ggtitle("")
b <- b + theme(legend.position="none") + ggtitle("")
c <- c + theme(legend.position="none") + ggtitle("")
d <- d + theme(legend.position="none") + ggtitle("")
e <- e + theme(legend.position="none") + ggtitle("")
f <- f + theme(legend.position="none") + ggtitle("")
g <- g + theme(legend.position="none") + ggtitle("")
h <- h + theme(legend.position="none") + ggtitle("")

pdf("pdf/subtyping/pcaWithoutLegends.pdf", height=5, width=12)
grid.arrange(a, b,c,d,e,f,g,h,nrow=2)
dev.off()

################################################################################
############################### UMAP VIZ #######################################
################################################################################
a <- ggscatter(merged_annotated, "UMAP1" ,"UMAP2", color="AD_status")+ scale_color_gradient(low="blue", high="red") + 
  ggtitle("Donor AD Status")
subgroup_cols <- c(
  "PathNonAD-TxMid" = "#fcbf49",
  "PathAD-TxAD" = "#ad2524",
  "PathAD-TxMid" = "#f38b20",
  "PathAD-TxNonAD" ="#90a955",
  "PathNonAD-TxNonAD" = "#264653",
  "PathNonAD-TxAD" = "#9d4edd"
)
b <- ggscatter(merged_annotated, "UMAP1", "UMAP2", color="subgroup", pal=subgroup_cols) + ggtitle("CPP-based Tx Subgroup")
cluster_cols <- c(
  "K1" = "#264653",
  "K2" = "#f38b20",
  "K3" = "#ad2524"
)
c <- ggscatter(merged_annotated,"UMAP1", "UMAP2", color="cluster_labels_for_big_heatmap", pal=cluster_cols) + ggtitle("CPP-based Tx Cluster")
d <- ggscatter(merged_annotated, "UMAP1", "UMAP2", color="tangles") + scale_color_gradient(low="gray", high="orange")
e <- ggscatter(merged_annotated, "UMAP1", "UMAP2", color="gpath") + scale_color_gradient(low="gray", high="darkgreen")
f <- ggscatter(merged_annotated, "UMAP1", "UMAP2", color="cogn_decline") + scale_color_gradient(low="gray", high = "purple")
g <- ggscatter(merged_annotated, "UMAP1", "UMAP2", color="age_death") + scale_color_gradient(low="gray", high = "brown")
h <- ggscatter(merged_annotated, "UMAP1", "UMAP2", color="msex") + scale_color_gradient(low="#4895ef", high = "#fb6f92")
grid.arrange(a,b,c,d,e,f,g,h,nrow=3)

################################################################################
############################## Correlations ####################################
################################################################################
metaVars <- c(
  "AD_status",
  "cogn_decline",
  "age_death",
  "msex",
  "amyloid",
  "plaq_n",
  "tangles",
  "gpath"
)

pcs_to_test <- paste0("PC", 1:10)
cor_matrix <- matrix(0, nrow = length(pcs_to_test), ncol = length(metaVars))
pval_matrix <- matrix(1, nrow = length(pcs_to_test), ncol = length(metaVars))
rownames(cor_matrix) <- pcs_to_test
colnames(cor_matrix) <- metaVars
rownames(pval_matrix) <- pcs_to_test
colnames(pval_matrix) <- metaVars
for (pc in pcs_to_test) {
  for (meta in metaVars) {
    x <- merged_annotated[[pc]]
    y <- merged_annotated[[meta]]
    test <- cor.test(x, y, method = "pearson", use = "complete.obs")
    cor_matrix[pc, meta] <- test$estimate
    pval_matrix[pc, meta] <- test$p.value
  }
}
annotation_matrix <- ifelse(pval_matrix < 0.001, "***",
                            ifelse(pval_matrix < 0.01, "**",
                                   ifelse(pval_matrix < 0.05, "*", "")))

# Apply bonf correction across all tests
padj_matrix <- matrix(
  p.adjust(as.vector(pval_matrix), method = "bonferroni"),
  nrow = nrow(pval_matrix),
  ncol = ncol(pval_matrix)
)
rownames(padj_matrix) <- rownames(pval_matrix)
colnames(padj_matrix) <- colnames(pval_matrix)

# Create annotation matrix based on adjusted p-values
annotation_matrix <- ifelse(padj_matrix < 0.001, "***",
                            ifelse(padj_matrix < 0.01, "**",
                                   ifelse(padj_matrix < 0.05, "*", "")))

pdf("pdf/subtyping/metaVarCor.pdf", height=4, width=7)
pheatmap(t(cor_matrix), 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = t(annotation_matrix),
         breaks = seq(-1, 1, length.out = 101),
         main = "PC vs MetaVar Correlations", 
         border_color = NA,
         number_color = "black", 
         fontsize_number=20
)
dev.off()

################################################################################
cellTypes <- get_fixed_celltype_order_for_plots()
endophenotypes <- c("AD_status", "amyloid", "tangles", "APOE_status", "plaq_n")
pcs_to_test <- paste0("PC", 1:5)

# First pass: collect all p-values for global BH correction
all_pvals <- c()
all_cors <- c()
test_info <- list()  # to track which test each p-value corresponds to

for (pheno in endophenotypes) {
  cellVars <- paste0("CPP.", pheno, ".", cellTypes)
  
  for (pc in pcs_to_test) {
    for (cell in cellVars) {
      x <- merged_annotated[[pc]]
      y <- merged_annotated[[cell]]
      test <- cor.test(x, y, method = "pearson", use = "complete.obs")
      
      all_pvals <- c(all_pvals, test$p.value)
      all_cors <- c(all_cors, test$estimate)
      test_info[[length(all_pvals)]] <- list(pheno = pheno, pc = pc, cell = cell)
    }
  }
}

# Apply bonf correction across all tests
all_padj <- p.adjust(all_pvals, method = "bonferroni")

# Second pass: build matrices with corrected p-values
heatmap_list <- list()

for (pheno in endophenotypes) {
  cellVars <- paste0("CPP.", pheno, ".", cellTypes)
  
  # initialize correlation and adjusted p-value matrices
  cor_matrix <- matrix(0, nrow = length(pcs_to_test), ncol = length(cellVars))
  padj_matrix <- matrix(1, nrow = length(pcs_to_test), ncol = length(cellVars))
  rownames(cor_matrix) <- pcs_to_test
  colnames(cor_matrix) <- cellVars
  rownames(padj_matrix) <- pcs_to_test
  colnames(padj_matrix) <- cellVars
  
  # fill matrices from stored results
  idx <- 1
  for (i in seq_along(test_info)) {
    info <- test_info[[i]]
    if (info$pheno == pheno) {
      cor_matrix[info$pc, info$cell] <- all_cors[i]
      padj_matrix[info$pc, info$cell] <- all_padj[i]
    }
  }
  
  # annotation stars based on adjusted p-values
  annotation_matrix <- ifelse(padj_matrix < 0.001, "***",
                              ifelse(padj_matrix < 0.01, "**",
                                     ifelse(padj_matrix < 0.05, "*", "")))
  
  # create the heatmap
  hm <- pheatmap(
    t(cor_matrix), 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_rownames = FALSE,
    legend = FALSE,
    border_color = NA,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    display_numbers = t(annotation_matrix),
    breaks = seq(-1, 1, length.out = 101),
    main = paste0(pheno),
    fontsize_number = 12
  )
  
  # store heatmap in list
  heatmap_list[[pheno]] <- hm
}

pdf("pdf/subtyping/cellType.PC.cor.pdf", height=7, width=10)
grid.arrange(heatmap_list[["AD_status"]]$gtable,
             heatmap_list[["amyloid"]]$gtable,
             heatmap_list[["tangles"]]$gtable,
             heatmap_list[["APOE_status"]]$gtable,
             heatmap_list[["plaq_n"]]$gtable, nrow=1)
dev.off()

# toggled show_rownames to TRUE for a moment to make this panel 
pdf("pdf/subtyping/cellType.PC.cor.withCellTypeLegend.pdf", height=7, width=6)
heatmap_list[["plaq_n"]]
dev.off()

################################################################################
######################### Plotting Cell Scores #################################
################################################################################

a <- ggscatter(merged_annotated, "PC1", "PC2", color="CPP.AD_status.Inh.ALCAM.TRPM3") + scale_color_gradient2(
  low = "blue", mid = "gray", high = "red",
  midpoint = 0, limits = c(-0.5, 0.5)
) + theme(legend.position="none") + ggtitle("Inh.ALCAM.TRPM3 CPP Score (AD Status)")
b <- ggscatter(merged_annotated, "PC1", "PC2", color="CPP.AD_status.Inh.L1.2.PAX6.SCGN") + scale_color_gradient2(
  low = "blue", mid = "gray", high = "red",
  midpoint = 0, limits = c(-0.5, 0.5)
) + theme(legend.position="none") + ggtitle("Inh.L1.2.PAX6.SCGN CPP Score (AD Status)")
c <- ggscatter(merged_annotated, "PC1", "PC2", color="CPP.AD_status.Inh.L1.PAX6.CA4") + scale_color_gradient2(
  low = "blue", mid = "gray", high = "red",
  midpoint = 0, limits = c(-0.5, 0.5)
) + theme(legend.position="none") + ggtitle("Inh.L1.PAX6.CA4 CPP Score (AD Status)")
d <- ggscatter(merged_annotated, "PC1", "PC2", color="CPP.AD_status.Oli") + scale_color_gradient2(
  low = "blue", mid = "gray", high = "red",
  midpoint = 0, limits = c(-0.5, 0.5)
) + theme(legend.position="none") + ggtitle("Oli CPP Score (AD Status)")
e <- ggscatter(merged_annotated, "PC1", "PC2", color="CPP.AD_status.Ast.DPP10") + scale_color_gradient2(
  low = "blue", mid = "gray", high = "red",
  midpoint = 0, limits = c(-0.5, 0.5)
) + theme(legend.position="none") + ggtitle("Ast.DPP10 CPP Score (AD Status)")
f <- ggscatter(merged_annotated, "PC1", "PC2", color="CPP.AD_status.Mic.P2RY12") + scale_color_gradient2(
  low = "blue", mid = "gray", high = "red",
  midpoint = 0, limits = c(-0.5, 0.5)
) + theme(legend.position="none") + ggtitle("MIc.P2RY12 CPP Score (AD Status)")


pdf("pdf/subtyping/cellScores.pdf", height=8, width=6)
grid.arrange(a,d,b,e,c,f,nrow=3)
dev.off()





my_breaks <- seq(-2, 2, length.out = 100)
my_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(my_breaks) - 1)
pheatmap(t(scale(combined_mat)),
         breaks = my_breaks,
         color = my_colors,
         fontsize_row = 4,
         show_colnames = FALSE,
         )
