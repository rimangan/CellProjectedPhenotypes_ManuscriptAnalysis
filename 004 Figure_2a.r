library(ggplot2)
library(ggridges)
library(ggrastr)
library(eimpute)

source("000 utils.r")

allcombined <- load_cell_scores()
meta <- load_meta_data()
## load cell type revised AD status
celltypeAll <- load_celltype_revised_phenotype("AD_status")


allcombinedsub <- allcombined[allcombined$projid %in% rownames(celltypeAll[["AD_status"]]),]
metasub <- meta[meta$projid %in% rownames(celltypeAll[["AD_status"]]),]
celltypeToShow <- "Mic.P2RY12"
allcombinedsub <- allcombinedsub[allcombinedsub$celltype == celltypeToShow,]
allcombinedsub$cellID <- rownames(allcombinedsub)
allcombinedsub$projid <- as.factor(allcombinedsub$projid)
allcombinedsub$cellID <- factor(allcombinedsub$cellID,levels = allcombinedsub[order(allcombinedsub$AD_status),"cellID"])
dfOrderProjid <- celltypeAll[["AD_status"]][,celltypeToShow]
names(dfOrderProjid) <- rownames(celltypeAll[["AD_status"]])
allcombinedsub$projid <- factor(allcombinedsub$projid,levels = names(sort(dfOrderProjid)))
dim(allcombinedsub)

## Figure 2a (donor-by-cell matrix)
p <- ggplot(allcombinedsub,aes(cellID,projid, fill = AD_status)) + geom_tile(color = "black", size = 0.25) +
            scale_fill_gradient2(low = "black", high = "black", mid = "black", midpoint = 0, limit = c(-2.1,2.1)) + theme_void() +
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none") #+

##with color
#p <- ggplot(allcombinedsub,aes(cellID,projid, col = AD_status)) + geom_point(size = 0.01) +
#            scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-2.1,2.1)) + theme_void() +
#            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),
#                  axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none") +
#            theme(panel.background = element_rect(fill = '#e8e8e8', colour = '#e8e8e8'))


allcombinedsub$diagnosis_individual <- metasub[match(allcombinedsub$projid,metasub$projid),"AD_status"] |> as.factor()
allcombinedsub <- allcombinedsub[!is.na(allcombinedsub$diagnosis_individual),]
allcombinedsub[allcombinedsub$diagnosis_individual == 1,"AD_status"] |> summary()
allcombinedsub[allcombinedsub$diagnosis_individual == 0,"AD_status"] |> summary()

## Figure 2a (distributions of cell scores of cells coming from non-AD and AD donors)
pdens <- ggplot(allcombinedsub,aes(x = AD_status,y = diagnosis_individual,fill = stat(x))) +
        geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2) + theme_minimal() +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-2,2)) +
         theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none") 
pdens <- rasterize(pdens,dpi = 200)

dfOrderProjid_df <- data.frame("projid" = names(dfOrderProjid),"AD_score" = unname(dfOrderProjid))
dfOrderProjid_df$projid <- factor(dfOrderProjid_df$projid,levels = names(sort(dfOrderProjid)))
dfOrderProjid_df$diagnosis <- meta[match(dfOrderProjid_df$projid,meta$projid),"AD_status"]
#dfOrderProjid_df$niareagansc <- meta[match(dfOrderProjid_df$projid,meta$projid),"niareagansc"]



## Figure 2a (cell type inferred AD phenotypes)
rightSide <- ggplot(dfOrderProjid_df,aes(x = "x",y = projid, fill = AD_score)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) + theme_minimal() + 
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")

## Figure 2a (donor-level phenotype)
leftside <- ggplot(dfOrderProjid_df,aes(x = "x",y = projid, fill = diagnosis)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5, limit = c(0,1)) + theme_minimal() + 
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")

## Figure 2a (distributions of celltype scores of non-AD and AD donors)
individualPdens <- ggplot(dfOrderProjid_df,aes(x = AD_score,y = as.factor(diagnosis),fill = stat(x))) +
        geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2) + theme_minimal() +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) + coord_flip()

ggsave(individualPdens, file = "Figures/Figure2/Fig2a_rightside_distributions.pdf",width = 4, height = 10)

top <- ggplot(allcombinedsub,aes(x = "x",y = cellID, fill = AD_status)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1.8,2)) + theme_minimal() + 
            theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none") + coord_flip()


head(allcombinedsub)
cellPdens <- ggplot(allcombinedsub,aes(x = AD_status,y = as.factor(diagnosis_individual),fill = stat(x))) +
        geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2) + theme_minimal() +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1.8,2))

ggsave(cellPdens, file = "Figures/Figure2/Fig2a_bottomside_distributions.pdf",width = 10, height = 2)

library(patchwork)
plot <- plot_spacer() + top + plot_spacer() +  leftside + p + rightSide + plot_layout(heights = c(1, 30), widths = c(1, 20, 1))
ggsave(plot, file = "Figures/Figure2/Fig2a_main.pdf",width = 20, height = 10)


