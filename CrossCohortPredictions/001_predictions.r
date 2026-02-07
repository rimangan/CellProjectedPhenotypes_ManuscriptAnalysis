library(ggplot2)
library(glmnet)

source("000_utils.r")
celltypeAll <- load_celltype_revised_phenotype(phenotypes = c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex"),
                                                dataset = "mathys")
predictedAll <- load_celltype_revised_phenotype(phenotypes = c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex"),
                                                      dataset = "jager")
meta <- load_BIG_meta_data()
phenotype <- "AD_status"
additional_vars <- c("msex", "pmi","age_death")# previously c("msex","APOE_status","pmi","age_death")
reference <- celltypeAll[[phenotype]] ## Our MIT data
target <- predictedAll[[phenotype]]   ## JAGER DATA
bothCells <- intersect(colnames(target),colnames(reference)) # ensure we analyze only cell types in both datasets
reference <- reference[,bothCells] # retaining cell types present in both sets
target <- target[,bothCells] # retaining cell types present in both sets
bothIDs <- intersect(rownames(reference),rownames(target))
## Using donors that were sampled in both cohorts to find features that are concordant
concordance <- sort(diag(cor(reference[bothIDs,],target[bothIDs,],method = "spearman")),decreasing = T)
concordanceplot <- data.frame(concordance)
concordanceplot$celltype <- rownames(concordanceplot)
concordanceplot$celltype <- factor(concordanceplot$celltype,levels = concordanceplot[order(concordanceplot$concordance),"celltype"])
concordanceplot$selected <- ifelse(concordanceplot$concordance>=0.5,"yes","no")
cp <- ggplot(concordanceplot,aes(concordance,celltype)) + geom_bar(stat = "identity") + theme_minimal() + labs(title = phenotype)

## Select cell types that are "concordant"
cellSelected <- names(concordance[concordance>=0.5])
reference <- reference[,cellSelected]
target <-   target[,cellSelected]
reference_matchingIDs <- reference[bothIDs,]
target_matchingIDs <- target[bothIDs,]

# diagonal of a correlation matrix between Jager and MIT shared donors based on CPP estimates,
# corresponding to correlation between the saem donor across datasets. We use IDs with positive correlations
# to learn a mapping between the MIT and Jager datasets
IDs_to_use <- diag(cor(t(reference_matchingIDs),t(target_matchingIDs)))
IDs_to_use <- names(IDs_to_use[IDs_to_use>0])

# use random forest wrapper to learn how to map reference CPP scores to target CPP scores
mapping_models <- learn_mapping(reference_matchingIDs[IDs_to_use,],target_matchingIDs[IDs_to_use,], ntree = 50, mtry = 15, nodesize = 10, maxnodes = 10,sampsize= length(IDs_to_use))

colnames(reference) <- make.names(colnames(reference))
colnames(target) <- make.names(colnames(target))
## apply mapping models on complete dataset
target <- align_data_sets(reference, target, mapping_models)


match_distributions <- function(reference,target){
    return(scale(target) * sd(reference) + mean(reference))
} 
target_matched <- lapply(1:ncol(reference),function(i){
    match_distributions(reference[,i],target[,i])
}) |> do.call(what = "cbind")

colnames(target_matched) <- colnames(target)
rownames(target_matched) <- rownames(target)

target <- target_matched |> as.data.frame()
target <- target[!rownames(target) %in% bothIDs,] |> as.data.frame()


#reference <- reference[!rownames(reference) %in% bothIDs,] |> as.data.frame()
reference$Target <- meta[match(rownames(reference),meta$projid),phenotype]
target$Target <- meta[match(rownames(target),meta$projid),phenotype]
target <- target[!is.na(target$Target),]
reference <- reference[!is.na(reference$Target),]
reference <- cbind(reference,meta[match(rownames(reference),meta$projid),c(additional_vars)])
target <- cbind(target,meta[match(rownames(target),meta$projid),c(additional_vars)])
target <- na.omit(target)
reference <- na.omit(reference)
target_IDs <- rownames(target)
reference_IDs <- rownames(reference)
combinedData <- rbind(reference,target)
baseModel <- additional_vars
fullModel <- colnames(combinedData)[colnames(combinedData) != "Target"]
CPP_Model <- fullModel[!fullModel %in% baseModel]

vars <- "baseModel"

# case for binary test variable (like AD_status)
perfomances <- lapply(c("baseModel","fullModel","CPP_Model"),function(model){   
    vars <- get(model) 
    Target_predicted <- lapply(target_IDs,function(x){
        message(x)
        trainData <- combinedData[!rownames(combinedData) %in% x,]|> as.data.frame()    
        test_data <- combinedData[!rownames(combinedData) %in% rownames(trainData),vars] |> as.data.frame()
        cv_model <- cv.glmnet(as.matrix(trainData)[,vars], trainData$Target, alpha = 0.5,family = "binomial")
        best_model <- glmnet(as.matrix(trainData)[,vars],trainData$Target,alpha = 0.5, lambda = cv_model$lambda.min, family = "binomial")
        return(predict(best_model,newx = as.matrix(test_data)[,vars]))
    }) |> do.call(what = "rbind") |> as.data.frame() 
    Target_predicted$true <- combinedData[target_IDs,"Target"]
    Reference_predicted <- lapply(reference_IDs,function(x){
        message(x)
        trainData <- combinedData[!rownames(combinedData) %in% x,]|> as.data.frame()    
        test_data <- combinedData[!rownames(combinedData) %in% rownames(trainData),vars] |> as.data.frame()
        cv_model <- cv.glmnet(as.matrix(trainData)[,vars], trainData$Target, alpha = 0.5,family = "binomial")
        best_model <- glmnet(as.matrix(trainData)[,vars],trainData$Target,alpha = 0.5, lambda = cv_model$lambda.min,family = "binomial")
        return(predict(best_model,newx = as.matrix(test_data)[,vars]))
    }) |> do.call(what = "rbind") |> as.data.frame() 
    Reference_predicted$true <- combinedData[reference_IDs,"Target"]
    Target_predicted$dataset <- "Jager"
    Target_predicted$projid <- target_IDs
    Reference_predicted$dataset <- "Kellis"
    Reference_predicted$projid <- reference_IDs
    toSave <- rbind(Target_predicted,Reference_predicted)
    saveRDS(toSave,sprintf("./Predictions/%s_%s.rds",phenotype,model))
    cor_target <- cor(Target_predicted$s0,Target_predicted$true,method = "spearman")
    cor_reference <- cor(Reference_predicted$s0,Reference_predicted$true,method = "spearman")
    Target_predicted$abs_errors <- abs(Target_predicted$s0 - Target_predicted$true)
    median_target <- median(Target_predicted$abs_errors)
    Reference_predicted$abs_errors <- abs(Reference_predicted$s0 - Reference_predicted$true)
    median_reference <- median(Reference_predicted$abs_errors)
    out <- c(cor_target,cor_reference,median_target,median_reference)
    return(out)    
})

# case for continuous variables (like amyloid, tangles, etc.)
source("000_utils.r")
celltypeAll <- load_celltype_revised_phenotype(phenotypes = c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex"),
                                               dataset = "mathys")
predictedAll <- load_celltype_revised_phenotype(phenotypes = c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex"),
                                                dataset = "jager")
meta <- load_BIG_meta_data()
phenotype <- "cogn_global_lv"
additional_vars <- c("msex", "pmi","age_death")# previously c("msex","APOE_status","pmi","age_death")
reference <- celltypeAll[[phenotype]] ## Our MIT data
target <- predictedAll[[phenotype]]   ## JAGER DATA
bothCells <- intersect(colnames(target),colnames(reference)) # ensure we analyze only cell types in both datasets
reference <- reference[,bothCells] # retaining cell types present in both sets
target <- target[,bothCells] # retaining cell types present in both sets
bothIDs <- intersect(rownames(reference),rownames(target))
## Using donors that were sampled in both cohorts to find features that are concordant
concordance <- sort(diag(cor(reference[bothIDs,],target[bothIDs,],method = "spearman")),decreasing = T)
concordanceplot <- data.frame(concordance)
concordanceplot$celltype <- rownames(concordanceplot)
concordanceplot$celltype <- factor(concordanceplot$celltype,levels = concordanceplot[order(concordanceplot$concordance),"celltype"])
concordanceplot$selected <- ifelse(concordanceplot$concordance>=0.5,"yes","no")
cp <- ggplot(concordanceplot,aes(concordance,celltype)) + geom_bar(stat = "identity") + theme_minimal() + labs(title = phenotype)

## Select cell types that are "concordant"
cellSelected <- names(concordance[concordance>=0.5])
reference <- reference[,cellSelected]
target <-   target[,cellSelected]
reference_matchingIDs <- reference[bothIDs,]
target_matchingIDs <- target[bothIDs,]

# diagonal of a correlation matrix between Jager and MIT shared donors based on CPP estimates,
# corresponding to correlation between the saem donor across datasets. We use IDs with positive correlations
# to learn a mapping between the MIT and Jager datasets
IDs_to_use <- diag(cor(t(reference_matchingIDs),t(target_matchingIDs)))
IDs_to_use <- names(IDs_to_use[IDs_to_use>0])

# use random forest wrapper to learn how to map reference CPP scores to target CPP scores
mapping_models <- learn_mapping(reference_matchingIDs[IDs_to_use,],target_matchingIDs[IDs_to_use,], ntree = 50, mtry = 15, nodesize = 10, maxnodes = 10,sampsize= length(IDs_to_use))

colnames(reference) <- make.names(colnames(reference))
colnames(target) <- make.names(colnames(target))
## apply mapping models on complete dataset
target <- align_data_sets(reference, target, mapping_models)


match_distributions <- function(reference,target){
  return(scale(target) * sd(reference) + mean(reference))
} 
target_matched <- lapply(1:ncol(reference),function(i){
  match_distributions(reference[,i],target[,i])
}) |> do.call(what = "cbind")

colnames(target_matched) <- colnames(target)
rownames(target_matched) <- rownames(target)

target <- target_matched |> as.data.frame()
target <- target[!rownames(target) %in% bothIDs,] |> as.data.frame()
#reference <- reference[!rownames(reference) %in% bothIDs,] |> as.data.frame()
reference$Target <- meta[match(rownames(reference),meta$projid),phenotype]
target$Target <- meta[match(rownames(target),meta$projid),phenotype]
target <- target[!is.na(target$Target),]
reference <- reference[!is.na(reference$Target),]
reference <- cbind(reference,meta[match(rownames(reference),meta$projid),c(additional_vars)])
target <- cbind(target,meta[match(rownames(target),meta$projid),c(additional_vars)])
target <- na.omit(target)
reference <- na.omit(reference)
target_IDs <- rownames(target)
reference_IDs <- rownames(reference)
combinedData <- rbind(reference,target)
baseModel <- additional_vars
fullModel <- colnames(combinedData)[colnames(combinedData) != "Target"]
CPP_Model <- fullModel[!fullModel %in% baseModel]

vars <- "baseModel"
perfomances <- lapply(c("baseModel","fullModel","CPP_Model"),function(model){   
  vars <- get(model) 
  Target_predicted <- lapply(target_IDs,function(x){
    message(x)
    trainData <- combinedData[!rownames(combinedData) %in% x,] |> as.data.frame()    
    test_data <- combinedData[!rownames(combinedData) %in% rownames(trainData), vars] |> as.data.frame()
    cv_model <- cv.glmnet(as.matrix(trainData)[,vars], trainData$Target, alpha = 0.5, family = "gaussian")
    best_model <- glmnet(as.matrix(trainData)[,vars], trainData$Target, alpha = 0.5, lambda = cv_model$lambda.min, family = "gaussian")
    return(predict(best_model, newx = as.matrix(test_data)[,vars]))
  }) |> do.call(what = "rbind") |> as.data.frame() 
  
  Target_predicted$true <- combinedData[target_IDs,"Target"]
  
  Reference_predicted <- lapply(reference_IDs,function(x){
    message(x)
    trainData <- combinedData[!rownames(combinedData) %in% x,] |> as.data.frame()    
    test_data <- combinedData[!rownames(combinedData) %in% rownames(trainData), vars] |> as.data.frame()
    cv_model <- cv.glmnet(as.matrix(trainData)[,vars], trainData$Target, alpha = 0.5, family = "gaussian")
    best_model <- glmnet(as.matrix(trainData)[,vars], trainData$Target, alpha = 0.5, lambda = cv_model$lambda.min, family = "gaussian")
    return(predict(best_model, newx = as.matrix(test_data)[,vars]))
  }) |> do.call(what = "rbind") |> as.data.frame() 
  
  Reference_predicted$true <- combinedData[reference_IDs,"Target"]
  
  Target_predicted$dataset <- "Jager"
  Target_predicted$projid <- target_IDs
  Reference_predicted$dataset <- "Kellis"
  Reference_predicted$projid <- reference_IDs
  
  toSave <- rbind(Target_predicted, Reference_predicted)
  saveRDS(toSave, sprintf("./Predictions/%s_%s.rds", phenotype, model))
  
  cor_target <- cor(Target_predicted$s0, Target_predicted$true, method = "spearman")
  cor_reference <- cor(Reference_predicted$s0, Reference_predicted$true, method = "spearman")
  Target_predicted$abs_errors <- abs(Target_predicted$s0 - Target_predicted$true)
  median_target <- median(Target_predicted$abs_errors)
  Reference_predicted$abs_errors <- abs(Reference_predicted$s0 - Reference_predicted$true)
  median_reference <- median(Reference_predicted$abs_errors)
  
  out <- c(cor_target, cor_reference, median_target, median_reference)
  return(out)    
})


library(ggpubr)
library(pROC)
library(gridExtra)
library(dplyr)

adStatusCppModel <- readRDS("Predictions/AD_status_CPP_Model.rds")
adStatusCppModel <- subset(adStatusCppModel, dataset == "Jager")
logit_model <- glm(true ~ s0, data = adStatusCppModel, family = binomial)
adStatusCppModel$predicted_prob <- predict(logit_model, type = "response")
t_test_result <- t.test(predicted_prob ~ true, data = adStatusCppModel)
a <- ggboxplot(adStatusCppModel, x = "true", y = "predicted_prob", 
               color = "true", fill = "true", alpha = 0.2, 
               palette = c("blue", "red"), add = "jitter", add.params=list(size=1)) +
  labs(x = "AD Status", y = "Predicted Probability of AD", title = "AD Status Prediction") +
  theme(legend.position = "none") +
  stat_compare_means(method = "t.test", label.x = 1.5)

### Amyloid ###
amyloidCppModel <- readRDS("Predictions/amyloid_CPP_Model.rds")
amyloidCppModel$Model <- "Inferred AD scores alone"
amyloidBaseModel <- readRDS("Predictions/amyloid_baseModel.rds")
amyloidBaseModel$Model <- "Metadata alone"
amyloidFullModel <- readRDS("Predictions/amyloid_fullModel.rds")
amyloidFullModel$Model <- "Inferred AD scores + Metadata"
mergedAmyloidData <- rbind(amyloidCppModel, amyloidBaseModel, amyloidFullModel)
mergedAmyloidData <- subset(mergedAmyloidData, dataset == "Jager")

amyloid_df_cor <- mergedAmyloidData %>%
  group_by(Model) %>%
  summarize(spearman_cor = cor(s0, true, method = "spearman")) 
amyloid_df_cor$Model <- factor(
  amyloid_df_cor$Model, 
  levels = c("Inferred AD scores alone",
             "Inferred AD scores + Metadata",
             "Metadata alone")
)
x <- ggbarplot(amyloid_df_cor, "Model", "spearman_cor", fill="Model", palette = c("#d62828", "#FFFFFF", "#0077b6"), label=TRUE, labs.nb.digits = 3)+labs(title="Amyloid Prediction", y="Spearman Correlation")+
  theme(axis.title.x=element_blank(), legend.position="none")

### Tangles ###
tanglesCppModel <- readRDS("Predictions/tangles_CPP_Model.rds")
tanglesCppModel$Model <- "Inferred AD scores alone"
tanglesBaseModel <- readRDS("Predictions/tangles_baseModel.rds")
tanglesBaseModel$Model <- "Metadata alone"
tanglesFullModel <- readRDS("Predictions/tangles_fullModel.rds")
tanglesFullModel$Model <- "Inferred AD scores + Metadata"
mergedTanglesData <- rbind(tanglesCppModel, tanglesBaseModel, tanglesFullModel)
mergedTanglesData <- subset(mergedTanglesData, dataset == "Jager")
tangles_df_cor <- mergedTanglesData %>%
  group_by(Model) %>%
  summarize(spearman_cor = cor(s0, true, method = "spearman"))
tangles_df_cor$Model <- factor(
  tangles_df_cor$Model, 
  levels = c("Inferred AD scores alone",
             "Inferred AD scores + Metadata",
             "Metadata alone")
)
y <- ggbarplot(tangles_df_cor, x = "Model", y = "spearman_cor", 
               fill = "Model", palette = c("#d62828", "#FFFFFF", "#0077b6"), 
               label = TRUE, labs.nb.digits = 3) +
  labs(title = "Tangle Density Prediction", 
       y = "Spearman Correlation") +
  theme(
    axis.title.x    = element_blank(),
    legend.position = "none"
  )

### Cognitive Impairment ###
cognCppModel <- readRDS("Predictions/cogn_global_lv_CPP_Model.rds")
cognCppModel$Model <- "Inferred AD scores alone"
cognBaseModel <- readRDS("Predictions/cogn_global_lv_baseModel.rds")
cognBaseModel$Model <- "Metadata alone"
cognFullModel <- readRDS("Predictions/cogn_global_lv_fullModel.rds")
cognFullModel$Model <- "Inferred AD scores + Metadata"
mergedCognData <- rbind(cognCppModel, cognBaseModel, cognFullModel)
mergedCognData <- subset(mergedCognData, dataset == "Jager")
cogn_df_cor <- mergedCognData %>%
  group_by(Model) %>%
  summarize(spearman_cor = cor(s0, true, method = "spearman"))
cogn_df_cor$Model <- factor(
  cogn_df_cor$Model, 
  levels = c("Inferred AD scores alone",
             "Inferred AD scores + Metadata",
             "Metadata alone")
)
z <- ggbarplot(cogn_df_cor, x = "Model", y = "spearman_cor", 
               fill = "Model", palette = c("#d62828", "#FFFFFF", "#0077b6"), 
               label = TRUE, labs.nb.digits = 3) +
  labs(title = "Cognitive Impairment Prediction", 
       y = "Spearman Correlation") +
  theme(
    axis.title.x    = element_blank(),
    legend.position = "none"
  )

pdf("threeModels.pdf", height=12, width=4)
grid.arrange(a,x,y,z, ncol=1)
dev.off()
