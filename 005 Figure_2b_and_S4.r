library(reshape2)
library(ggplot2)
library(patchwork)

source("000 utils.r")

## Load meta data
meta <- load_meta_data()
## load cell type revised phenotypes
celltypeAll <- load_celltype_revised_phenotype(c("AD_status","amyloid","tangles","CR","APOE_status","age_death","plaq_n","cogn_global_lv","msex"))
### get fixed celltype order ###
celltypesorder <- get_fixed_celltype_order_for_plots()

##Flip 1= male to female = 1
#celltypeAll[["msex"]] <- celltypeAll[["msex"]] *-1
#celltypeAll[["cogn_global_lv"]] <- celltypeAll[["cogn_global_lv"]] *-1

## Figure 2b (cell type importance)
pheno_to_show <- "cogn_global_lv"
combined <- get_revised_combined(pheno_to_show,bin = F)
MFI <- cor(celltypeAll[[pheno_to_show]],combined[,1],method ="spearman",use = "pairwise.complete.obs") |> melt()
MFI$Var1 <- factor(MFI$Var1,levels = celltypesorder)
importancePlot <- ggplot(MFI,aes(Var1,range01(value))) + geom_bar(stat = "identity") + theme_minimal() +
                            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none") +
                            labs(y = 'Cell type importance')
## Figure 2b (cell type revised phenotypes heatmap)
phenoInferred <- celltypeAll[[pheno_to_show]]
phenoInferred <- melt(as.matrix(phenoInferred))
phenoInferred$Var1 <- as.character(phenoInferred$Var1)
phenoInferred$Var1 <- factor(phenoInferred$Var1,levels = names(combined[order(combined[,1]),]))
phenoInferred$Var2 <- factor(phenoInferred$Var2,levels = celltypesorder)


limit <- c(min(phenoInferred$value),max(phenoInferred$value))
ranges <- 0.3
limit <-  c(-ranges,ranges)
phenoInferred[phenoInferred$value < -ranges,"value"] <- -ranges
phenoInferred[phenoInferred$value > ranges,"value"] <- ranges

mainP <- ggplot(phenoInferred,aes(Var2,Var1, fill = value)) + geom_tile() +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = limit) + theme_minimal() + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
            

## Figure 2b (donor-level phenotype)
meta_copy <- meta
meta_copy$projid<- factor(meta_copy$projid,levels = names(combined[order(combined[,1]),]))
lim <- getLimforPlot(scale(meta_copy[,pheno_to_show]))
leftside <- ggplot(meta_copy,aes(x = "x",y = projid, fill = scale(get(pheno_to_show)))) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-lim,lim)) + theme_minimal() + 
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")

## Figure 2b (combined inferred phenotype)
meta_copy <- meta_copy[match(rownames(combined),meta_copy$projid),]
meta_copy$combinedScore <- scale(combined[,1])
lim <- getLimforPlot(meta_copy$combinedScore)
rightside <- ggplot(meta_copy,aes(x = "x",y = projid, fill = combinedScore)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-lim,lim)) + theme_minimal() + 
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")

## combine 
combined_plots <- plot_spacer() + importancePlot + plot_spacer() + leftside + mainP + rightside +  plot_layout(widths = c(1, 52, 1),heights = c(2,5))
ggsave(combined_plots, file = "Figures/Figure2/Fig2b/Cognitive_impairment.pdf",width = 10, height = 5)



meta <- meta[match(rownames(celltypeAll[["amyloid"]]),meta$projid),]
names(celltypeAll)
out <- cor(celltypeAll[["amyloid"]],get_revised_combined("amyloid",F,celltypes = "all"),method = "spearman") |> as.data.frame()

out <- out[order(out[,1]),,drop = F]



############ SUPP 4 ###############
######## From Judy ################
library(gridExtra)
library(dplyr)
cppScores <- readRDS(file = "Data/all_cell_scores_new3.rds")
source("000 utils.r")
## Load meta data
meta <- load_meta_data()
meta$AD_status_donor <- ifelse(meta$niareagansc %in% c(4,3),0,1) ## Binarize AD status
meta$msex_donor <- meta$msex
meta <- meta %>% mutate(projid = as.integer(projid))

# defining major cell types
Exc <- cppScores$celltype[grepl("^Exc\\.", cppScores$celltype)] |> unique()
Inh <- cppScores$celltype[grepl("^Inh\\.", cppScores$celltype)] |> unique()
Ast <- cppScores$celltype[grepl("^Ast\\.", cppScores$celltype)] |> unique()
Mic <- c("Mic.P2RY12", "Mic.MKI67")
ImmuneVasc <- c("T.cells", "CAMs", "Per", "SMC", "End", "Fib.FLRT2", "Fib.SLC4A4")
Oli <- c("Oli")
OPC <- c("OPCs")

# Create the major cell types list
majorCellTypes <- list(
  Exc = Exc,
  Inh = Inh,
  Ast = Ast,
  Mic = Mic,
  ImmuneVasc = ImmuneVasc,
  Oli = Oli,
  OPC = OPC
)

cppScores <- cppScores %>% 
  left_join(meta %>% select(projid, AD_status_donor, msex_donor), by = "projid")

calculate_donor_averages <- function(cell_types, major_type_name, cpp_score_type) {
  # Filter data for the specific cell types
  filtered_data <- cppScores[cppScores$celltype %in% cell_types, ]
  
  if(nrow(filtered_data) == 0) {
    warning(paste("No data found for", major_type_name))
    return(NULL)
  }
  
  # Calculate average CPP score per donor for the specified CPP score type
  donor_averages <- filtered_data %>%
    group_by(projid) %>%
    summarise(
      mean_cpp_score = mean(.data[[cpp_score_type]], na.rm = TRUE),
      n_cells = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      major_cell_type = major_type_name,
      cpp_score_type = cpp_score_type
    )
  
  return(donor_averages)
}

all_donor_averages <- list()

# Define the CPP score types to analyze (based on Figure 2 of the paper)
cpp_score_types <- c("AD_status", "amyloid", "tangles", "cogn_global_lv", "plaq_n", "APOE_status", "CR", "msex")
cpp_score_labels <- c("AD Status", "Amyloid", "Tangle Density", "Cognitive Impairment", "Neuritic Plaque", "APOE Status", "Cognitive Resilience", "Biological Sex (M/F)")

counter <- 1
for(i in 1:length(majorCellTypes)) {
  major_type_name <- names(majorCellTypes)[i]
  cell_types <- majorCellTypes[[i]]
  
  print(paste("Processing", major_type_name, "with", length(cell_types), "cell subtypes"))
  
  # Calculate for each CPP score type
  for(j in 1:length(cpp_score_types)) {
    cpp_score_type <- cpp_score_types[j]
    cpp_label <- cpp_score_labels[j]
    
    donor_avg <- calculate_donor_averages(cell_types, major_type_name, cpp_score_type)
    if(!is.null(donor_avg)) {
      all_donor_averages[[counter]] <- donor_avg
      counter <- counter + 1
    }
  }
}

# cell type order
lvl_major <- c("Exc","Inh","Ast","Oli","OPC","Mic","ImmuneVasc")

# Combine all results
combined_averages <- bind_rows(all_donor_averages)
final_data <- combined_averages %>%
  left_join(
    cppScores %>% select(projid, AD_status_donor, msex_donor) %>% distinct(),
    by = "projid"
  )
final_data <- final_data %>%
  mutate(AD_status_label = case_when(
    AD_status_donor == 0 ~ "Non-AD",
    AD_status_donor == 1 ~ "AD", 
    TRUE ~ NA_character_
  )) %>%
  mutate(msex_label = case_when(
    msex_donor == 0 ~ "Female",
    msex_donor == 1 ~ "Male",
    TRUE ~ NA_character_
  )) %>%
  mutate(major_cell_type = factor(major_cell_type, 
                                  levels = rev(lvl_major)))


pAD <- ggboxplot(subset(final_data, cpp_score_type == "AD_status"), x="major_cell_type", y="mean_cpp_score", color="AD_status_label", fill="AD_status_label", alpha=0.5, add="jitter") +
  scale_color_manual(values= c("AD" = "#ec2426", "Non-AD" = "#3b54a4")) + geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="none")+labs(title="AD Status")
pAmyloid <- ggboxplot(subset(final_data, cpp_score_type == "amyloid"), x="major_cell_type", y="mean_cpp_score", color="AD_status_label", fill="AD_status_label", alpha=0.5, add="jitter") +
  scale_color_manual(values= c("AD" = "#ec2426", "Non-AD" = "#3b54a4")) + geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="none")+labs(title="Amyloid")
pTangles <- ggboxplot(subset(final_data, cpp_score_type == "tangles"), x="major_cell_type", y="mean_cpp_score", color="AD_status_label", fill="AD_status_label", alpha=0.5, add="jitter") +
  scale_color_manual(values= c("AD" = "#ec2426", "Non-AD" = "#3b54a4")) + geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="none")+labs(title="Tangle Density")
pCognDecline <- ggboxplot(subset(final_data, cpp_score_type == "cogn_global_lv"), x="major_cell_type", y="mean_cpp_score", color="AD_status_label", fill="AD_status_label", alpha=0.5, add="jitter") +
  scale_color_manual(values= c("AD" = "#ec2426", "Non-AD" = "#3b54a4")) + geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="none")+labs(title="Cognitive Impairment")
pPlaque <- ggboxplot(subset(final_data, cpp_score_type == "plaq_n"), x="major_cell_type", y="mean_cpp_score", color="AD_status_label", fill="AD_status_label", alpha=0.5, add="jitter") +
  scale_color_manual(values= c("AD" = "#ec2426", "Non-AD" = "#3b54a4")) + geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="none")+labs(title="Neuritic Plaque")
pAPOE <- ggboxplot(subset(final_data, cpp_score_type == "APOE_status"), x="major_cell_type", y="mean_cpp_score", color="AD_status_label", fill="AD_status_label", alpha=0.5, add="jitter") +
  scale_color_manual(values= c("AD" = "#ec2426", "Non-AD" = "#3b54a4")) + geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="none")+labs(title="APOE Status")
pCR <- ggboxplot(subset(final_data, cpp_score_type == "CR"), x="major_cell_type", y="mean_cpp_score", color="AD_status_label", fill="AD_status_label", alpha=0.5, add="jitter") +
  scale_color_manual(values= c("AD" = "#ec2426", "Non-AD" = "#3b54a4")) + geom_hline(yintercept=0, linetype="dashed") + coord_flip() + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="none")+labs(title="Cognitive Reslience")
pdf("Figures/Supplementary/Supp4_Mean_Boxplots/S4a.pdf", height=7, width=24)
grid.arrange(pAD, pAmyloid, pTangles, pCognDecline, pPlaque, pCR, nrow=1)
dev.off()



# Sanity check: make S4a with cell type labels shown

pAD_y <- ggboxplot(subset(final_data, cpp_score_type == "AD_status"), 
                   x="major_cell_type", y="mean_cpp_score", 
                   color="AD_status_label", fill="AD_status_label", 
                   alpha=0.5, add="jitter") +
  scale_color_manual(values=c("AD"="#ec2426","Non-AD"="#3b54a4")) +
  geom_hline(yintercept=0, linetype="dashed") +
  coord_flip() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(),   # show y title
        axis.text.y = element_text(),  # show y labels
        legend.position="none") +
  labs(title="AD Status", y="Cell type")

pAmyloid_y <- pAD_y %+% subset(final_data, cpp_score_type == "amyloid") + labs(title="Amyloid")
pTangles_y <- pAD_y %+% subset(final_data, cpp_score_type == "tangles") + labs(title="Tangle Density")
pCognDecline_y <- pAD_y %+% subset(final_data, cpp_score_type == "cogn_global_lv") + labs(title="Cognitive Impairment")
pPlaque_y <- pAD_y %+% subset(final_data, cpp_score_type == "plaq_n") + labs(title="Neuritic Plaque")
pAPOE_y <- pAD_y %+% subset(final_data, cpp_score_type == "APOE_status") + labs(title="APOE Status")
pCR_y <- pAD_y %+% subset(final_data, cpp_score_type == "CR") + labs(title="Cognitive Resilience")

pdf("Figures/Supplementary/Supp4_Mean_Boxplots/S4_with_cell_type_y_axis.pdf", height=7, width=28)
grid.arrange(pAD_y, pAmyloid_y, pTangles_y, pCognDecline_y, pPlaque_y, pAPOE_y, pCR_y, nrow=1)
dev.off()

 