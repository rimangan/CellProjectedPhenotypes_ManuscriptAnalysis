load_celltype_revised_phenotype <- function(phenotypes, dataset = "mathys"){
    require(eimpute)
    require(reshape2)
    allcombined <- load_cell_scores(dataset)
    celltypeAll <- lapply(phenotypes,function(y){
        message(y)
        celltypes <- unique(allcombined$celltype) |> as.character()
        cors <- lapply(celltypes,function(x){
            message(x)
            sub <- allcombined[allcombined$celltype == x,]
            sub_pc_ADCells <- lapply(unique(sub$projid),function(ID){        
                mean(sub[sub$projid == ID,y],na.rm=T) 
            }) |> do.call(what = "rbind")
            sub_pc_ADCells <- as.data.frame(sub_pc_ADCells)
            sub_pc_ADCells$projid <- unique(sub$projid)
            return(sub_pc_ADCells)
        })
        allIDs <- lapply(cors,function(cor){
            return(cor$projid)
        }) |> unlist() |> unique()
        all_AD_PCs <- lapply(cors,function(cor){
            cor[match(allIDs,cor$projid),1]
        }) |> do.call(what = "cbind")
        rownames(all_AD_PCs) <- allIDs
        colnames(all_AD_PCs) <- celltypes    
        all_AD_PCs <- all_AD_PCs[,colSums(apply(all_AD_PCs,2,function(x)!is.na(x))) >=150]
        all_AD_PCs <- all_AD_PCs[colSums(apply(all_AD_PCs,1,function(x)!is.na(x))) >= 15,]    
        imputed <- eimpute(all_AD_PCs, 5)
        imputed <- imputed[[1]]
        rownames(imputed) <- rownames(all_AD_PCs)
        colnames(imputed) <- colnames(all_AD_PCs)
        all_AD_PCs <- imputed |> as.data.frame()
        return(all_AD_PCs)
    })
    names(celltypeAll) <- phenotypes
    return(celltypeAll)
}

load_cell_scores <- function(dataset = "mathys"){
    if(dataset == "mathys"){
        ## Load cell scores
        allcombined <- readRDS(file = "Data/all_cell_scores_new3.rds")
        #allcombined <- readRDS(file = "Data/all_cell_scores.rds")
        ## Assign column names
        columnNames <- c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv",
                    "amyloid_Res","tangles_Res","CR_Res","APOE_multi_Res","APOE_status_Res","age_death_Res","plaq_n_Res","cogn_global_lv_Res",
                    "msex","msex_Res")
       # colnames(allcombined) <- c(columnNames,"projid","celltypes")
        #colnames(allcombined) <- c("AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex","projid","celltypes")
        out <- allcombined
    }else if(dataset == "jager"){
        allJager <- readRDS(file = "Data/all_cell_scores_jager_new.rds")        
        #colnames(allJager) <- c("cellBarcode","projid","celltypes","AD_status","amyloid","tangles","CR","APOE_multi","APOE_status","age_death","plaq_n","cogn_global_lv","msex")
        out <- allJager
    }
    return(out)    
}

load_BIG_meta_data <- function(){
    require(openxlsx)
    metaall <- read.csv("Data/dataset_1369_basic.csv")
    metaall$AD_status <- ifelse(metaall$niareagansc %in% c(4,3),0,1)
    metaall$APOE_status <- ifelse(metaall$apoe_genotype %in% c(34,44),1,0)
    metaall$projid <-as.character(metaall$projid)
    metaall$cogn_global_lv <- metaall$cogn_global_lv  *-1
    metaall$msex <- (metaall$msex - 1) *-1
    metaall$APOE_multi <- ifelse(metaall$apoe_genotype == 44,2,NA)    
    metaall$APOE_multi <- ifelse(metaall$apoe_genotype == 34,1,metaall$APOE_multi) 
    metaall$APOE_multi <- ifelse(metaall$apoe_genotype == 33,0,metaall$APOE_multi)
    metaall$APOE_multi <- ifelse(metaall$apoe_genotype == 23,-1,metaall$APOE_multi)
    metaall$APOE_multi <- ifelse(metaall$apoe_genotype == 22,-2,metaall$APOE_multi)
    metaall$APOE_multi <- ifelse(metaall$apoe_genotype == 24,0.5,metaall$APOE_multi)    
    metaall[is.na(metaall$APOE_multi),"APOE_multi"] <- 0
    metaall$plaq_n <- log(metaall$plaq_n+1)
    metaall$amyloid <- log(metaall$amyloid+1)
    metaall$tangles <- log(metaall$tangles+1)
    return(metaall)
}

load_meta_data <- function(){
    ##Load meta data##
    meta <- read.csv("Data/Metadata.csv")
    meta$AD_status <- ifelse(meta$niareagansc %in% c(4,3),0,1)
    ## Calculate Cognitive resilience 
    meta$CR <- lm(cogn_global_lv~gpath,data = meta)$residuals
    ## Calculate Cognitive Decline Resilience 
    CDRs <- lm(cogng_random_slope~gpath,data = meta)$residuals
    meta[names(CDRs),"CDR"] <- CDRs
    meta$APOE_status <- ifelse(meta$apoe_genotype %in% c(34,44),1,0) |> as.numeric()
    meta$APOE_multi <- ifelse(meta$apoe_genotype == 44,2,NA)    
    meta$APOE_multi <- ifelse(meta$apoe_genotype == 34,1,meta$APOE_multi) 
    meta$APOE_multi <- ifelse(meta$apoe_genotype == 33,0,meta$APOE_multi)
    meta$APOE_multi <- ifelse(meta$apoe_genotype == 23,-1,meta$APOE_multi)
    meta$APOE_multi <- ifelse(meta$apoe_genotype == 22,-2,meta$APOE_multi)
    meta$APOE_multi <- ifelse(meta$apoe_genotype == 24,0.5,meta$APOE_multi)
    meta$plaq_n <- log(meta$plaq_n+1)
    meta$amyloid <- log(meta$amyloid+1)
    meta$tangles <- log(meta$tangles+1)
    meta$cogn_global_lv <- meta$cogn_global_lv  *-1
    return(meta)
}

get_fixed_celltype_order_for_plots <- function(){
    celltypeAll_specific <- load_celltype_revised_phenotype("AD_status")[[1]]    
    ID_order <- names(sort(rowSums(celltypeAll_specific)))
    celltypeAll_specific_long <- melt(as.matrix(celltypeAll_specific))
    celltypes <- unique(celltypeAll_specific_long$Var2)
    Inh <- celltypes[grepl("Inh",celltypes)]|> as.character()
    Exc <- celltypes[grepl("Exc",celltypes)]|> as.character()
    Ast <- celltypes[grepl("Ast",celltypes)]|> as.character()
    Mic <- celltypes[grepl("Mic",celltypes)]|> as.character()
    Oli <- celltypes[grepl("Oli",celltypes)]|> as.character()
    OPC <- celltypes[grepl("OPC",celltypes)]|> as.character()
    Tcells <- celltypes[grepl("T.cells",celltypes)]|> as.character()
    Fib <- celltypes[grepl("Fib.FLRT2",celltypes)]|> as.character()
    End <- celltypes[grepl("End",celltypes)]|> as.character()
    Per <- celltypes[grepl("Per",celltypes)]|> as.character()
    SMC <- celltypes[grepl("SMC",celltypes)]|> as.character()
    CAM <- celltypes[grepl("CAMs",celltypes)]|> as.character()
    celltypesorder <- c(names(sort(colSums(celltypeAll_specific[,Exc]),decreasing = T)),
                        names(sort(colSums(celltypeAll_specific[,Inh]),decreasing = T)),
                        names(sort(colSums(celltypeAll_specific[,Ast]),decreasing = T)),
                        Oli,
                        OPC,
                        names(sort(colSums(celltypeAll_specific[,Mic]),decreasing = T)),
                        Tcells,CAM,Per,SMC,End,Fib)|> as.character()
    return(celltypesorder)
}

get_revised_combined <- function(pheno,bin,celltypes = "all"){
    require(glmnet)
    if(celltypes[1] == "all"){
        cellinferred <- celltypeAll[[pheno]]
    }else{
        cellinferred <- celltypeAll[[pheno]][,celltypes]
    }    
    meta <- meta[match(rownames(cellinferred),meta$projid),]
    if(bin){
        cv_model <- cv.glmnet(as.matrix(cellinferred), meta[,pheno], alpha = 1, family = "binomial")
        best_model <- glmnet(as.matrix(cellinferred), meta[,pheno],alpha = 1, lambda = cv_model$lambda.min, family = "binomial")
    }else{
        if(pheno %in% c("amyloid","plaq_n","tangles")){
            cv_model <- cv.glmnet(as.matrix(cellinferred), log(meta[,pheno]+1), alpha = 0)
            best_model <- glmnet(as.matrix(cellinferred), log(meta[,pheno]+1),alpha = 0, lambda = cv_model$lambda.min)
        }else{
            cv_model <- cv.glmnet(as.matrix(cellinferred), meta[,pheno], alpha = 0)
            best_model <- glmnet(as.matrix(cellinferred), meta[,pheno],alpha = 0, lambda = cv_model$lambda.min)
        }        
    }    
    probs <- predict(best_model, type = "response", newx = as.matrix(cellinferred))
    return(probs)
}
range01 <- function(x){
    out <- (x-median(x))/(max(x)-median(x))
    out[out<=0] <-0
    return(out)
}

ceiling_dec <- function(x, level=1){
    round(x + 5*10^(-level-1), level)
} 
getLimforPlot <- function(vec){
    lim <- max(abs(c(min(vec),max(vec)))) |> ceiling_dec()
    return(lim)
}

prepareData_progressionPlots <- function(path, celltypes,IDs){
    out <- lapply(celltypes,function(x){
        plotData <- data.frame(meta[IDs,path],AD_components[IDs,x],x)
        colnames(plotData) <- c("var1","var2","var3")
        plotData <- plotData[order(plotData$var1),]
        plotData$index <- 1:nrow(plotData)
        return(plotData)
    }) |> do.call(what = "rbind")
    return(out)
}


train_validation_test_split <- function(data){
    n <- nrow(data)
    # Create indices for splitting
    train_indices <- sample(1:n, size = 0.7 * n) # 70% for training
    remaining_indices <- setdiff(1:n, train_indices)
    validation_indices <- sample(remaining_indices, size = 0.5 * length(remaining_indices)) # 15% for validation
    test_indices <- setdiff(remaining_indices, validation_indices) # Remaining 15% for testing
    # Subset the data
    train_data <- data[train_indices, ]
    validation_data <- data[validation_indices, ]
    test_data <- data[test_indices, ]
    return(list(train = train_data,validation = validation_data , test = test_data))
}

perform_grid_search <- function(train_data, validation_data, grid){
    for (i in 1:nrow(grid)) {    
        # Extract parameters for this iteration
        ntree <- grid$ntree[i]
        mtry <- grid$mtry[i]
        nodesize <- grid$nodesize[i]
        maxnodes <- grid$maxnodes[i]
        sampsize <- grid$sampsize[i]
        replace <- grid$replace[i]
        n_PCs <- grid$PCs[i]
        # Error handling using tryCatch
        tryCatch({
            # Train the model
            rf_model <- randomForest(Target ~ ., data = train_data[,c(paste0("PC",1:n_PCs),"Target")], 
                                    #Target ~ .,
                                    #data = train_data,
                                    ntree = ntree, 
                                    mtry = mtry, 
                                    nodesize = nodesize, 
                                    maxnodes = ifelse(is.na(maxnodes), NULL, maxnodes), 
                                    sampsize = sampsize, 
                                    replace = replace)
            
            # Predict on validation set
            validation_predictions <- predict(rf_model, newdata = validation_data[,c(paste0("PC",1:n_PCs),"Target")])
            #validation_predictions <- predict(rf_model, newdata = validation_data)
            # Compute RMSE on validation set
            validation_actuals <- validation_data$Target
            validation_rmse <- sqrt(mean((validation_predictions - validation_actuals)^2))
            # Store the RMSE
            results$Validation_RMSE[i] <- validation_rmse
            results$Validation_cor[i] <- cor(validation_predictions,validation_actuals,method = "spearman")                            
            #results$Jager_cor[i] <- cor(validation_data[grepl("Jager",rownames(validation_data)),"Target"], validation_predictions[grepl("Jager",names(validation_predictions))],method = "spearman")
            #results$MIT_cor[i] <- cor(validation_data[grepl("MIT",rownames(validation_data)),"Target"], validation_predictions[grepl("MIT",names(validation_predictions))],method = "spearman")
            results$Jager_median_abs_error[i] <- median(abs(validation_data[grepl("Jager",rownames(validation_data)),"Target"] - validation_predictions[grepl("Jager",names(validation_predictions))]))
            results$MIT_median_abs_error[i] <- median(abs(validation_data[grepl("MIT",rownames(validation_data)),"Target"] - validation_predictions[grepl("MIT",names(validation_predictions))]))
            #message(results$Jager_cor[i])
        }, error = function(e) {
            # Handle error: Print message and skip this iteration
            cat("Error at iteration", i, "with parameters:", 
                "ntree =", ntree, 
                "mtry =", mtry, 
                "nodesize =", nodesize, 
                "maxnodes =", maxnodes, 
                "sampsize =", sampsize, 
                "replace =", replace, "\n")
            cat("Error message:", e$message, "\n")
            # Leave Validation_RMSE as NA for this iteration
        })
    }    
    grid <- results[order(grid$Validation_cor,decreasing = TRUE),]
    return(grid)
}

learn_mapping <- function(reference, target, ntree, mtry, nodesize, maxnodes,sampsize){
    require(randomForest)
    rownames(reference) <- paste0("ref_",rownames(reference))
    rownames(target) <- paste0("target_",rownames(target))
    R <- reference
    T <- target
    models <- vector("list", ncol(R))
    for (j in 1:ncol(R)) {
        # Prepare data frame for randomForest
        df <- data.frame(R_col = R[,j], T)
        # Fit a random forest model predicting R[,j] from all columns of T
        models[[j]] <- randomForest(R_col ~ ., data=df, ntree=ntree,mtry = mtry,nodesize=nodesize, maxnodes=maxnodes, sampsize=sampsize)
    }
    names(models) <- colnames(reference)
    return(models)
}

align_data_sets <- function(reference,target, models){
    newtarget <- lapply(models,function(model){predict(model,target)}) |> do.call(what = "cbind")
    colnames(newtarget) <- colnames(target)
    target <- as.data.frame(newtarget)
    return(target)
}

euclidean_distance <- function(row1, row2) {
  sqrt(sum((row1 - row2)^2))
}

ExpandingNeighborhoodLocalSearch <- function(eucs, allcombined, phenotype, start_cell, n_start = 5, n_max = 200, temperature = 0.1){
  storing <- c(start_cell)  
  repeat {
    #message("Current cell: ", start_cell)
    
    own_score <- allcombined[start_cell, phenotype]
    
    found_next <- FALSE  # track if we found any valid neighbor in [n_start..n_max]
    
    # Attempt from n_start up to n_max
    for (N in n_start:n_max) {
      # 1) Find the N nearest neighbors by distance
      nearest_n <- sort(eucs[start_cell, ]) |> head(n = N)
      
      # 2) Exclude already visited
      nearest_n <- nearest_n[!names(nearest_n) %in% storing]
      
      # 3) Keep only neighbors with strictly lower (or slightly lower) score
      neighbor_ok <- allcombined[names(nearest_n), phenotype] < (own_score)
      nearest_n   <- nearest_n[neighbor_ok]
      
      # If we have valid neighbors, pick one
      if (length(nearest_n) > 0) {
        # 4) Compute delta scores
        delta_scores <- own_score - allcombined[names(nearest_n), phenotype]
        names(delta_scores) <- names(nearest_n)
        
        # 5) Weight them by softmax(distance * -1) * delta_scores
        weights <- softmax(nearest_n * -1, temperature = temperature) * delta_scores      

        # 6) Pick the neighbor with the best (largest) weight
        pick_sorted <- sort(weights)
        next_cell   <- tail(pick_sorted, 1)
        next_cell   <- names(next_cell)  # get the name, not the numeric value
        
        # 7) If the chosen neighbor is the same cell, break 
        if (next_cell == start_cell) {
          message("Next cell would be the same as current; stopping.")
          found_next <- FALSE  # effectively force a stop
        } else {
          # We found a valid next cell
          found_next <- TRUE
        }
        
        break  # break out of the for(N in ...) loop
      }
    }
    
    # If we didn't find any valid neighbor up to n_max, we're done
    if (!found_next) {
      message("No more strictly lower-scoring neighbors (up to N = ", n_max, "). Stopping.")
      break
    }
    
    # Otherwise, append this next cell to the path and reset N to n_start
    storing    <- c(storing, next_cell)
    start_cell <- next_cell
    
    # The loop repeats with N reset back to n_start
  }

  # 'storing' now contains your path
  return(storing)
}


fix_celltype_name <- function(celltype) {
  celltype <- switch(celltype,
    "End" = "Endothelial",
    "OPC" = "OPCs",
    "Oli" = "Oligodendrocytes",
    celltype  # Default case: return the original celltype
  )
  return(celltype)
}


load_metabolite_hmdb_ids_mapping <- function(){
    files <- list.files("/tudelft.net/staff-umbrella/brainQTLs/Gerard/AD430/Everything/Data")
    files <- files[grepl("primary_secondary",files)]

    meta_ids <- lapply(files,function(file){
        readRDS(sprintf("/tudelft.net/staff-umbrella/brainQTLs/Gerard/AD430/Everything/Data/%s",file))
    }) |> do.call(what = "rbind")


    all_meta_ids <- apply(meta_ids,1,function(row){
        primary <- row[[2]]
        secondaries <- row[[3]]
        secondaries <- strsplit(secondaries,split = "; ") |> unlist()
        data.frame("primary" = primary,"secondaries" = c(primary,secondaries))
    }) |> do.call(what = "rbind")
    all_meta_ids <- na.omit(all_meta_ids)
    return(all_meta_ids)
}

load_variant_data <- function(){
    vars <- readRDS("/tudelft.net/staff-umbrella/brainQTLs/Gerard/AD430/ROSMAP_SEAAD_AD_107.rds")
    varData <- vars[[2]]
    varInfo <- vars[[1]]
    samplesheet <- read.csv2("/tudelft.net/staff-umbrella/brainQTLs/Gerard/AD430/Samplesheet.csv")
    projID_individualID_mapping <- samplesheet[,c("projid","individualID")]
    #projID_individualID_mapping <- projID_individualID_mapping[projID_individualID_mapping$projid %in% meta$projid,]
    projID_individualID_mapping <- projID_individualID_mapping[projID_individualID_mapping$individualID %in% colnames(varData),]
    varData <- varData[,projID_individualID_mapping$individualID]
    colnames(varData) <- projID_individualID_mapping$projid
    return(varData)
}
