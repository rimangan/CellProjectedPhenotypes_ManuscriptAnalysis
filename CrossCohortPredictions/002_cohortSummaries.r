library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(circlize)


source("000_utils.r")
celltypeAll <- load_celltype_revised_phenotype(phenotypes = c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex"),
                                               dataset = "mathys")
predictedAll <- load_celltype_revised_phenotype(phenotypes = c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex"),
                                                dataset = "jager")

tsaiKellis_adStatus <- as.data.frame(celltypeAll["AD_status"])
deJager_adStatus <- as.data.frame(predictedAll["AD_status"])

tsaiKellis_adStatus$projid <- rownames(tsaiKellis_adStatus)
deJager_adStatus$projid <- rownames(deJager_adStatus)

meta <- load_BIG_meta_data()

joined_df <- full_join(
  tsaiKellis_adStatus %>% mutate(in_TsaiKellis = TRUE),
  deJager_adStatus %>% mutate(in_deJager = TRUE),
  by = "projid"
) %>%
  mutate(
    in_TsaiKellis = replace_na(in_TsaiKellis, FALSE),
    in_deJager = replace_na(in_deJager, FALSE),
    Cohort = case_when(
      in_TsaiKellis & in_deJager ~ "Both",
      in_TsaiKellis & !in_deJager ~ "TsaiKellisOnly",
      !in_TsaiKellis & in_deJager ~ "DeJagerOnly"
    )
  ) %>%
  select(projid, Cohort)

meta_unique <- meta %>% distinct(projid, .keep_all = TRUE)
merged_df <- left_join(joined_df, meta_unique, by = "projid")

metabolites <- readRDS("../Data/metabolites_fully_processed.rds")
has_metabolomics <- metabolites$projid
merged_df <- merged_df %>%
  mutate(HasMetabolomics = projid %in% has_metabolomics)


################################################################################
################ Heatmap #######################################################
################################################################################
library(ComplexHeatmap)
library(RColorBrewer)

# Define the annotations
annotations <- merged_df %>%
  select(Cohort, HasMetabolomics, AD_status, gpath, msex, age_death)

cohort_key <- factor(
  annotations$Cohort,
  levels = c("TsaiKellisOnly", "Both", "DeJagerOnly")
)

metab_key <- with(annotations, ifelse(
  Cohort == "TsaiKellisOnly",
  ifelse(HasMetabolomics, 2, 1),   # FALSE first, TRUE second
  ifelse(
    Cohort == "Both",
    ifelse(HasMetabolomics, 1, 2), # TRUE first, FALSE second
    0                              # DeJagerOnly: neutral
  )
))
gpath_key <- annotations$gpath
reordered_indices <- order(
  cohort_key,
  metab_key,
  gpath_key
)

#reordered_indices <- order(annotations$Cohort, annotations$HasMetabolomics, annotations$gpath)
annotations <- annotations[reordered_indices, ]

# Create a list to store the annotations
annotation_list <- list()

# Add 'Cohort' annotation
cohort_colors <- c("Both" = "#ee9b00", "DeJagerOnly" = "#0a9396", "TsaiKellisOnly" = "#ae2012")
annotation_list[["Cohort"]] <- HeatmapAnnotation(
  Cohort = annotations$Cohort,
  col = list(Cohort = cohort_colors),
  show_legend = TRUE
)

# Add 'gpath' annotation
gpath_color_func <- colorRamp2(c(min(annotations$gpath), max(annotations$gpath)), c("#fefae0", "#bc6c25"))
annotation_list[["gpath"]] <- HeatmapAnnotation(
  gpath = annotations$gpath,
  col = list(gpath = gpath_color_func),
  show_legend = TRUE
)

# Add 'AD_status' annotation
ad_status_colors <- c("0" = "blue", "1" = "red")
annotation_list[["AD_status"]] <- HeatmapAnnotation(
  AD_status = annotations$AD_status,
  col = list(AD_status = ad_status_colors),
  show_legend = TRUE
)


# Add 'msex' annotation
msex_colors <- setNames(c("#219ebc", "lightpink"), unique(annotations$msex))
annotation_list[["msex"]] <- HeatmapAnnotation(
  msex = annotations$msex,
  col = list(msex = msex_colors),
  show_legend = TRUE
)


# Add 'gpath' annotation
age_death_color_func <- colorRamp2(c(min(annotations$age_death), max(annotations$age_death)), c("#fefae0", "#283618"))
annotation_list[["age_death"]] <- HeatmapAnnotation(
  age_death = annotations$age_death,
  col = list(age_death = age_death_color_func),
  show_legend = TRUE
)

# Add 'HasMetabolomics' annotation
has_metabolomics_colors <- c("TRUE" = "black", "FALSE" = "lightgray")
annotation_list[["HasMetabolomics"]] <- HeatmapAnnotation(
  HasMetabolomics = annotations$HasMetabolomics,
  col = list(HasMetabolomics = has_metabolomics_colors),
  show_legend = TRUE
)

ha <- HeatmapAnnotation(
  Cohort = annotations$Cohort,
  gpath = annotations$gpath,
  AD_status = annotations$AD_status,
  msex = annotations$msex,
  age_death = annotations$age_death,
  HasMetabolomics = annotations$HasMetabolomics,
  col = list(
    Cohort = cohort_colors,
    gpath = gpath_color_func,
    HasMetabolomics = has_metabolomics_colors,
    AD_status = ad_status_colors,
    msex = msex_colors,
    age_death = age_death_color_func
  ),
  show_legend = TRUE,
  annotation_legend_param = list(
    Cohort = list(nrow = 3),
    gpath = list(title = "gpath"),
    age_death = list(title = "age_death")
  )
)

matrix_data <- matrix(1, nrow = 1, ncol = nrow(annotations))
ht <- Heatmap(
  matrix_data,
  top_annotation = ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = c("white"),
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

pdf("pdf/donorHeatmap.pdf", height=10, width=5)
draw(ht)
dev.off()
