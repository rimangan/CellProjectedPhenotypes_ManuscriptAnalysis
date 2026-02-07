library(data.table)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(gridExtra)

#setwd("/tudelft.net/staff-umbrella/brainQTLs/Gerard/AD430/Everything/Riley") #specific to gerard's local

#Load meta data
meta <- read.csv("Data/Metadata.csv")
## Binarize AD status
meta$AD_status <- ifelse(meta$niareagansc %in% c(4,3),0,1)
## Calculate Cognitive resilience 
meta$CR <- lm(cogn_global_lv~gpath,data = meta)$residuals
##Load cell scores
#allcombined <- readRDS(file = "Data/all_cell_scores.rds")
allcombined <- readRDS(file = "Data/all_cell_scores_new3.rds")

## Assign column names
#columnNames <- c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv",
#               "amyloid_Res","tangles_Res","CR_Res","APOE_multi_Res","APOE_status_Res","age_death_Res","plaq_n_Res","cogn_global_lv_Res",
#               "msex","msex_Res")

#columnNames <- c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex")           
#colnames(allcombined) <- c(columnNames,"projid","celltypes")

## Load UMAP coordinates
BigUMAP <- fread("Data/UMAP_coordinates.coords.tsv")

## Get barcodes that match barcodes from all combined
splitted <- strsplit(BigUMAP$V1,"[_]")
splitted <- sapply(splitted,function(x)x[3])
BigUMAP$barcode <- splitted

### Clean barcodes
cleaned_strings <- gsub("-\\d+$", "", rownames(allcombined))
allcombined$barcode <- cleaned_strings
BigUMAP <- as.data.frame(BigUMAP)
allcombined$UMAP_1 <- BigUMAP[match(allcombined$barcode,BigUMAP$barcode),"V2"]
allcombined$UMAP_2 <- BigUMAP[match(allcombined$barcode,BigUMAP$barcode),"V3"]
allcombined <- allcombined[!is.na(allcombined$UMAP_1),]
allcombined <- na.omit(allcombined)

## Get UMAP coordinates for labels
annoData <- lapply(unique(allcombined$celltype),function(h){
        colMeans(allcombined[allcombined$celltype == h,c("UMAP_1","UMAP_2")])
    }) |> do.call(what = "rbind") |> as.data.frame()
annoData$label <- unique(allcombined$celltype)


## Cap AD scores for visualization purposes
allcombined$AD_status_capped <- allcombined$AD_status
allcombined[allcombined$AD_status_capped>=0.35,"AD_status_capped"] <- 0.35
allcombined[allcombined$AD_status_capped<=-0.35,"AD_status_capped"] <- -0.35

## Add original donor level AD diagnosis
allcombined$AD_diagnosis <- meta[match(allcombined$projid,meta$projid),"AD_status"]

## Random down sample for visualization purposes
celltypes <- unique(allcombined$celltype)
downSampled <- lapply(celltypes,function(cell){
    message(cell)
    sub <- allcombined[allcombined$celltype == cell,]
    nTotal <- nrow(sub)
    n_to_sample <- round(nTotal/10) * 0.95
    return(sub[sample(rownames(sub),n_to_sample),])
}) |> do.call(what ='rbind')


allcombined$celltype <- factor(allcombined$celltype,levels = sample(unique(allcombined$celltype)))

## Figure 1a
Pcols <- ggplot(allcombined,aes(UMAP_1,UMAP_2, col = celltype)) + geom_point(size = 0.002,alpha = 0.1)  + theme_void() +
 theme(legend.position = "none") #+ geom_text_repel(data = annoData, aes(x = UMAP_1, y = UMAP_2, label = label),size = 5, col = "black")
## Rasterize
#Pcols <- rasterize(Pcols,dpi = 100)
ggsave(Pcols, file = "Fig1a_without_labels.pdf",width = 20, height = 20)

## Figure 1b
Pcols <- ggplot(allcombined,aes(UMAP_1,UMAP_2,col = as.factor(AD_diagnosis))) + theme_void()  +geom_point(size = 0.002,alpha = 1) +
    scale_color_manual(values=c("#0059FF","#FF1818"))  + theme(legend.position = "none")
## Rasterize
#Pcols <- rasterize(Pcols,dpi = 100)
ggsave(Pcols, file = "Figures/Figure1/Fig1b.pdf",width = 20, height = 20)

## Figure 1c
Pcols <- ggplot(data = allcombined,aes(UMAP_1,UMAP_2, col = AD_status_capped)) + geom_point(size = 0.002,alpha = 1) +
        scale_color_gradient2(low = "#0059FF", high = "#FF1818", mid = "#E8E1DA",
        midpoint = 0, limit = c(-0.35,0.35), space = "Lab",name="phenotype") +
        theme_void() + theme(legend.position = "top") + 
        theme(legend.position = "none")
#Pcols <- rasterize(Pcols,dpi = 100)
ggsave(Pcols, file = "Figures/Figure1/Fig1c.pdf",width = 20, height = 20)

## Rasterize
Pcols <- rasterize(Pcols,dpi = 100)



#### SUPP FIGURE 2 #####
plot_umap_feature <- function(df, feature, midpoint = 0, cap = 1) {
  df <- df %>%
    mutate(feature_value = pmax(pmin(!!sym(feature), cap), -cap))  # cap values between -1 and 1
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = feature_value)) +
    geom_point(size=0.05, alpha=1) +
    scale_color_gradient2(
      low = "#0059FF", high = "#FF1818", mid = "#E8E1DA",
      midpoint = midpoint, limits = c(-cap, cap)) +
    theme_classic() +
    theme(
      legend.position = "none",
    )
}
features <- c("AD_status", "amyloid", "tangles", "CR", "APOE_status", "plaq_n", "cogn_global_lv", "msex")

for (feat in features) {
  currPlot <- plot_umap_feature(allcombined, feat)
  png(
    filename = paste0("Figures/Supplementary/", feat, ".png"),
    width = 4, height = 4, units = "in", res = 300
  )
  print(currPlot)
  dev.off()
}


#### SUPP 3B (From Gulfem)####
library(dplyr)
library(pheatmap)

cppScores <- readRDS("Data/all_cell_scores_new3.rds")
cppScores <- cppScores %>% mutate(projid = as.integer(as.character(projid)))
cppEndo <- cppScores %>%
  select(
    AD_status = AD_status,
    amyloid = amyloid,
    APOE_status = APOE_status,
    cognitive_resilience = CR,
    cognitive_impairment = cogn_global_lv,
    sex = msex,
    neuritic_plaque = plaq_n,
    tangle_density = tangles,
    celltype = celltype,
  )

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

cellTypeList <- lapply(majorCellTypes, function(subtypes) {
  cppEndo[cppEndo$celltype %in% subtypes, ]
})

endophenotypes <- c("AD_status", "amyloid", "APOE_status", "cognitive_resilience", "cognitive_impairment", "sex",
                    "neuritic_plaque", "tangle_density")
nice_labels <- c(
  "AD_status" = "AD Status",
  "amyloid" = "Amyloid",
  "tangle_density" = "Tangle Density",
  "cognitive_resilience" = "Cognitive Resilience",
  "APOE_status" = "APOE Status",
  "neuritic_plaque" = "Neuritic Plaques",
  "cognitive_impairment" = "Cognitive Impairment",
  "sex" = "Sex"
)

plot_correlation_heatmap <- function(df, title, cluster = TRUE) {
  mat <- df[, endophenotypes]
  mat <- data.frame(lapply(mat, function(x) as.numeric(as.character(x))))
  cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  
  rownames(cor_mat) <- nice_labels[rownames(cor_mat)]
  colnames(cor_mat) <- nice_labels[colnames(cor_mat)]
  breaks <- seq(-1, 1, length.out = 101)
  pheatmap(cor_mat,
           cluster_rows = cluster,
           cluster_cols = cluster,
           display_numbers = FALSE,
           main = title,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = breaks,
           fontsize_number = 10,
           fontsize = 12,
           border_color = NA,
           number_color = "black")
}

pEnh <- plot_correlation_heatmap(cellTypeList[["Exc"]], "Excitatory Neurons")
pInh <- plot_correlation_heatmap(cellTypeList[["Inh"]], "Inhibitory Neurons")
pAst <- plot_correlation_heatmap(cellTypeList[["Ast"]], "Astrocytes")
pMic <- plot_correlation_heatmap(cellTypeList[["Mic"]], "Microglia")
pImmuneVasc <- plot_correlation_heatmap(cellTypeList[["ImmuneVasc"]], "Immune & Vascular")
pOli <- plot_correlation_heatmap(cellTypeList[["Oli"]], "Oligodendrocytes")
pOPC <- plot_correlation_heatmap(cellTypeList[["OPC"]], "OPCs")

pdf("Figures/Supplementary/Supp2_CorrUmap/corrPlot.pdf", height=24, width=14)
grid.arrange(pEnh$gtable, pInh$gtable, pAst$gtable, pMic$gtable, pImmuneVasc$gtable, pOli$gtable, pOPC$gtable, ncol=2)
dev.off()
