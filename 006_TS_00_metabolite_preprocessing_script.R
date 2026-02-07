# ----------------------------- Packages ---------------------------------------------------------------------
library(eimpute)
library(ggplot2)
library(stringr)
library(patchwork)
library(ggrepel)

# ----------------------------- Utils & metadata -------------------------------------------------------------
source("000 utils.r")
meta <- load_BIG_meta_data()   
cpp_scores <- readRDS("Data/combined_metadata.rds")

# ensure projid uniqueness:
dup_n <- sum(duplicated(cpp_scores$projid))
if (dup_n > 0) {
  message(sprintf("[cpp_scores] Found %d duplicate projid(s); keeping first per projid.", dup_n))
  cpp_scores <- cpp_scores[!duplicated(cpp_scores$projid), , drop = FALSE]
}
rownames(cpp_scores) <- cpp_scores$projid

# ----------------------------- Inputs -----------------------------------------------------------------------
meta_dict   <- read.csv("Data/ROSMAP Metabolon HD4 Data Dictionary.csv")
# --- Drop all unnamed metabolites from meta_dict ---
meta_dict <- subset(meta_dict, TYPE != "UNNAMED")

refmet      <- read.csv("Data/refmet.csv")  
metabo_data <- read.csv("Data/ROSMAP_Metabolon_HD4_Brain514_assay_data.csv")
projID_individualID_mapping <- read.csv2("Data/Samplesheet.csv")
projID_individualID_mapping <- projID_individualID_mapping[, c("projid","individualID")]
pc2_groups_AD <- readRDS("Data/projid_pc2_groups_ADonly.rds")  
merged_df_final <-readRDS("Data/merged_df_final.rds")

# Map individualID -> projid
metabo_data$projid <- projID_individualID_mapping[
  match(metabo_data$individualID, projID_individualID_mapping$individualID),
  "projid"
]

# Keep all biological samples that have a valid projid
metabo_data <- metabo_data[!is.na(metabo_data$projid) & nzchar(metabo_data$projid), , drop = FALSE]

# remaining duplicate projids?
dup_n <- sum(duplicated(metabo_data$projid))
if (dup_n > 0) {
  message(sprintf("[metabo_data] Found %d duplicate projid(s); keeping first per projid.", dup_n))
  metabo_data <- metabo_data[!duplicated(metabo_data$projid), , drop = FALSE]
}

rownames(metabo_data) <- metabo_data$projid

# ----------------------------- Metabolite column selection --------------------------------------------------
# Identify "no-name" metabolites in dict
noName_metabolites <- meta_dict[grepl("^X\\d+$", meta_dict$SHORT_NAME), "CHEM_ID"]
noName_metabolites <- paste0("X", noName_metabolites)

# robust selection: columns that look like X<digits> and appear in dict
met_cols_regex   <- grepl("^X\\d+$", colnames(metabo_data))
met_cols_in_dict <- gsub("^X", "", colnames(metabo_data)) %in% meta_dict$CHEM_ID
met_cols <- which(met_cols_regex & met_cols_in_dict)

metabo_dataSel <- metabo_data[, met_cols, drop = FALSE]
metabo_dataSel <- metabo_dataSel[, !(colnames(metabo_dataSel) %in% noName_metabolites), drop = FALSE]

# ----------------------------- Missingness filter ----------------------------------
colKeep <- colSums(is.na(metabo_dataSel)) <= nrow(metabo_dataSel) * 0.25
metabo_dataSel <- metabo_dataSel[, colKeep, drop = FALSE]

# ----------------------------- Impute (raw scale) and PQN normalization, then log transform -----------------
set.seed(1)  
metabo_dataSel_imputed <- eimpute(metabo_dataSel, 5)[[1]]

# PQN normalization: compute reference spectrum as the median across samples
ref_sample <- apply(metabo_dataSel_imputed, 2, median, na.rm = TRUE)
metabo_dataSel_norm <- t(apply(metabo_dataSel_imputed, 1, function(sample) {
  ratio <- sample / ref_sample
  factor <- median(ratio, na.rm = TRUE)
  sample / factor
}))

# Log-transform normalized data
metabo_dataSel_im <- log(metabo_dataSel_norm)

# plot distribution change across these preprocessing steps
library(ggplot2)
v_raw      <- as.numeric(unlist(metabo_dataSel))                 
v_imputed  <- as.numeric(unlist(metabo_dataSel_imputed))         
v_pqn      <- as.numeric(unlist(metabo_dataSel_norm))            
v_logfinal <- as.numeric(unlist(metabo_dataSel_im))              

plot_raw      <- log1p(v_raw[!is.na(v_raw)])
plot_imputed  <- log1p(v_imputed)
plot_pqn      <- log1p(v_pqn)
plot_logfinal <- v_logfinal  

df_steps <- rbind(
  data.frame(val = plot_raw,      step = "Raw (log1p)"),
  data.frame(val = plot_imputed,  step = "Imputed (log1p)"),
  data.frame(val = plot_pqn,      step = "PQN (log1p)"),
  data.frame(val = plot_logfinal, step = "Final log")
)

ggplot(df_steps, aes(x = val, fill = step)) +
  geom_histogram(bins = 100, position = "identity", alpha = 0.35) +
  labs(x = "log-scale intensity (comparable units)", y = "Count",
       title = "Distributions across preprocessing steps") +
  theme_minimal(base_size = 13)


# Harmonize row/col names
rownames(metabo_dataSel_im) <- metabo_data$projid
colnames(metabo_dataSel_im) <- colnames(metabo_dataSel)

# ----------------------------- Create and save "metabolites_fully_processed" (projid + SHORT_NAME columns) =====

stopifnot(exists("metabo_dataSel_im"),
          is.matrix(metabo_dataSel_im) || is.data.frame(metabo_dataSel_im),
          exists("meta_dict"),
          all(c("CHEM_ID","SHORT_NAME") %in% names(meta_dict)))

# 1) Pull CHEM_IDs from current columns (strip leading "X")
chem_from_cols <- gsub("^X", "", colnames(metabo_dataSel_im))

# 2) Build a mapping: current col -> SHORT_NAME
map_df <- data.frame(
  col_old     = colnames(metabo_dataSel_im),
  CHEM_ID     = chem_from_cols,
  stringsAsFactors = FALSE
)
map_df$SHORT_NAME <- meta_dict$SHORT_NAME[match(map_df$CHEM_ID, meta_dict$CHEM_ID)]

# 3) Fallback: if a SHORT_NAME is missing/blank, keep CHEM_ID as the name
map_df$SHORT_NAME[is.na(map_df$SHORT_NAME) | !nzchar(map_df$SHORT_NAME)] <- map_df$CHEM_ID

# 4) Preserve true SHORT_NAMEs, only enforce uniqueness 
raw_labels <- map_df$SHORT_NAME
raw_labels <- trimws(raw_labels)

# Build final names ensuring the first entry "projid" is reserved, and any duplicates (including "projid" among metabolites) get ".1", ".2", ...
final_names <- make.unique(c("projid", raw_labels), sep = ".")

changed_idx <- which(final_names[-1] != raw_labels)
if (length(changed_idx)) {
  message("Deduplicated metabolite names (original -> unique):")
  for (i in head(changed_idx, 20)) {
    message(sprintf("  %s -> %s", raw_labels[i], final_names[i + 1]))
  }
  if (length(changed_idx) > 20) message(sprintf("  ... and %d more", length(changed_idx) - 20))
}

# 5) Construct final df 
metabolites_fully_processed <- as.data.frame(metabo_dataSel_im, check.names = FALSE)
metabolites_fully_processed$projid <- rownames(metabolites_fully_processed)
metabolites_fully_processed <- metabolites_fully_processed[, c(ncol(metabolites_fully_processed),
                                                               1:(ncol(metabolites_fully_processed)-1))]

# 6) Apply new column names
colnames(metabolites_fully_processed) <- final_names

# 7) sanity checks
stopifnot(!anyDuplicated(metabolites_fully_processed$projid))
stopifnot(!anyDuplicated(colnames(metabolites_fully_processed)[-1]))

# 8) save
saveRDS(metabolites_fully_processed, file = "Data/metabolites_fully_processed.rds")

# message
message(sprintf("[metabolites_fully_processed] %d samples × %d metabolites (+projid).",
                nrow(metabolites_fully_processed),
                ncol(metabolites_fully_processed) - 1))
