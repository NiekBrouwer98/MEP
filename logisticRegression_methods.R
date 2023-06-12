# This file contains all necessary methods to run the logistic regression experiments

predictorInitialFile <- function(){
  # Structures
  TMEStructures <- here('scratch/ptLeveLTMEInteractions.fst')
  TMEStructures <- read_fst(TMEStructures, as.data.table = T)
  IDs <- getCells()[, .(ImageNumber, isDiscovery)]
  IDs <- IDs[, .SD[1], by = ImageNumber]
  TMEStructures <- merge(x = TMEStructures, y = IDs, by = 'ImageNumber') 
  TMEStructures[, isTestCohort := ifelse(isDiscovery, 'train', 'test')]
  setnames(TMEStructures, 
           c('nInteractionsPerStructure', 'TotalInteractions'),
           c('module_count', 'total_communities'))
  TMEStructures[, proportion := module_count / total_communities]
  TMEStructures[, proportion := gtools::logit(ifelse(proportion == 1, proportion - 1e-6, proportion + 1e-6))]
  structures <- getStructures()
  predictors <- TMEStructures[, .(ImageNumber, TMEStructure, proportion, total_communities, isTestCohort)]
  predictors[, weights := sum(unique(total_communities)), by = ImageNumber]
  predictors[, total_communities := NULL]
  predictors <- dcast(predictors, ImageNumber + weights + isTestCohort ~ TMEStructure, value.var = 'proportion')
  oldNames <- grep('ImageNumber|weights|isTestCohort', names(predictors), invert = TRUE, value = TRUE)
  newNames <- paste0('TMEStructure', oldNames)
  setnames(predictors, oldNames, newNames)
  
  # Network properties
  nwP <- here('scratch/NetworkProperties.rds')
  nwP <- read_rds(nwP)
  standardiseNwP <- function(name){
    dt <- nwP[[name]]
    setnames(dt, grep('communities_', names(dt), value = TRUE), 'communityID')
    dt[, type := name]
    return(dt)
  }
  nwP <- rbindlist(lapply(c('stroma', 'tumour'), standardiseNwP))
  IDs <- getCells()[, .(ImageNumber, metabric_id)][, .SD[1], by = ImageNumber]
  nwP <- merge(x = nwP, y = IDs, by = 'ImageNumber') 
  nwP[, assortativity := NULL] #has NAs
  measureVars <- grep('ImageNumber|communityID|type|metabric_id', names(nwP), 
                      invert = TRUE, value = TRUE)
  PtmeasureVars <- paste0(measureVars, 'PtMean')
  nwP[, eval(PtmeasureVars) := lapply(.SD, mean), by = .(type, ImageNumber), .SDcols = measureVars]
  nwP <- nwP[, .SD[1], by = .(type, ImageNumber)]
  nwP <- nwP[, .SD, .SDcols = c('ImageNumber', 'type', PtmeasureVars)]
  nwP <- melt(nwP, id.vars = c('ImageNumber', 'type'))
  nwP[, variable := paste0(variable, '_', type)]
  nwP <- dcast(nwP, ImageNumber ~ variable)
  nwP <- na.omit(nwP)
  nwPredictorsTumour <- grep('_tumour', names(nwP), invert = F, value = T)
  nwPredictorsTME <- grep('_stroma', names(nwP), invert = F, value = T)
  predictors <- merge(x = predictors, y = nwP, by = 'ImageNumber') 
  
  cells <- getCellsAlternative()[(isTumour)]
  
  # Compute proportions: all, tumour, stroma, vascular, interface
  mkProportionBy <- function(byvars, dt, suffix){
    
    countsBy <- byvars
    totalsBy <- setdiff(byvars, 'meta_description')
    outVar <- paste0('proportion_',suffix)
    
    dt[, counts := .N, by = byvars]
    dt[, totals := .N, by = totalsBy]
    dt[, eval(outVar) := (counts / totals)]
    dt[, eval(c('counts', 'totals')) := NULL]
    return(dt)
  }
  
  ptCellVars <- c('ImageNumber', 'meta_description') 
  mkProportionBy(byvars = ptCellVars,dt = cells, suffix = 'all')
  
  mkProportionBy(
    byvars = c(ptCellVars, 'is_epithelial'),
    dt = cells, suffix = 'isEpi')
  
  # cells[, is_vascular := (Parent_vessel != 0)]
  # mkProportionBy(
  # 	byvars = c(ptCellVars, 'is_epithelial', 'is_vascular'),
  # 	dt = cells, suffix = 'isVesselByEpi')
  # 
  # mkProportionBy(
  # 	byvars = c(ptCellVars, 'is_epithelial', 'is_interface'),
  # 	dt = cells, suffix = 'isInterfaceByEpi')
  
  proportionVars <- grep('proportion_', names(cells), value = T)
  phenotypeVars <- grep('meta_|phenotype|colours', names(cells), value = T)
  indicatorVars <- grep('^is_', names(cells), value = T)
  indicatorVars <- setdiff(indicatorVars, c('is_normal', 'is_dcis'))
  
  toKeep <- c(ptCellVars, proportionVars, phenotypeVars, indicatorVars)
  proportionsOut <- cells[, .SD, .SDcols = toKeep]
  proportionsOut <- melt(proportionsOut,
                         id.vars = c('ImageNumber', phenotypeVars, indicatorVars),
                         measure.vars = proportionVars,
                         value = 'proportion',
                         variable = 'type'
  )	
  proportionsOut[, type := gsub('proportion_', '', type)]
  proportionsOut[type == 'isEpi', 
                 type := ifelse(is_epithelial, 'tumour', 'stroma')]
  # proportionsOut[type == 'isVesselByEpi', 
  # 	type := ifelse(is_vascular, 'vesselByEpi', 'NotVesselByEpi')]
  # proportionsOut[type == 'isInterfaceByEpi', 
  # 	type := ifelse(is_interface, 'interfaceByEpi', 'NotInterfaceByEpi')]
  # proportionsOut[grep('ByEpi',type), newSuffix := gsub('[0-9]*', '', meta_description)]
  # proportionsOut[, type := gsub('ByEpi', '', type)]
  # proportionsOut[!is.na(newSuffix), type := paste0(type, newSuffix)][,
  # 	newSuffix := NULL]
  proportionsOut <- proportionsOut[, .SD[1], by = .(ImageNumber, meta_description, type)]
  proportionsOut[, type := gsub(' ', '', type)]
  
  proportionsOut[, check := sum(proportion), by = .(ImageNumber, type)]
  stopifnot(all.equal(rep(1, nrow(proportionsOut)), proportionsOut[['check']]))
  proportionsOut[, check := NULL]
  
  outfile <- here('scratch/cellPhenotypeProportionsAlternative.fst')
  write_fst(proportionsOut, outfile)
  
  # Cell phenotype proportions
  cellPhenotypes <- here('scratch/cellPhenotypeProportionsAlternative.fst')
  cellPhenotypes <- read_fst(cellPhenotypes, as.data.table = T)
  cellPhenotypes <- cellPhenotypes[grep('tumour|stroma', type)]
  cellPhenotypes <- cellPhenotypes[, .(ImageNumber, type, meta_description, proportion)]
  allCombinations <- adt(expand.grid(unique(cellPhenotypes[,ImageNumber]), unique(cellPhenotypes[,meta_description])))
  setnames(allCombinations, c('ImageNumber', 'meta_description'))
  allCombinations <- merge(x = allCombinations, 
                           y = cellPhenotypes[!duplicated(meta_description), .(meta_description, type)], 
                           by = 'meta_description') 
  cellPhenotypes <- merge(x = allCombinations, y = cellPhenotypes, 
                          by = c('ImageNumber', 'meta_description', 'type'), all.x = T) 
  # cellPhenotypes[, meta_description := paste0(meta_description, '_CPh')]
  cPh_tumour <- cellPhenotypes[type == 'tumour', unique(meta_description)]
  cPh_tme <- cellPhenotypes[type == 'stroma', unique(meta_description)]
  cellPhenotypes[is.na(proportion), proportion := 0]
  cellPhenotypes[, checkBothCompartments := sum(proportion), by = .(ImageNumber, type)]
  cellPhenotypes <- cellPhenotypes[checkBothCompartments > 0][, checkBothCompartments := NULL]
  cellPhenotypes <- dcast(cellPhenotypes, ImageNumber ~ meta_description, value.var = 'proportion')
  cellPhenotypes <- na.omit(cellPhenotypes) # samples that contain both tumour and stromal cells
  predictors <- merge(x = predictors, y = cellPhenotypes, by = 'ImageNumber')
  saveRDS(predictors, file = here('scratch/predictors_intfile.rds'))
  
  # Add original cell type densities
  cells <- getCells()[(isTumour)]
  ptCellVars <- c('ImageNumber', 'meta_description') 
  mkProportionBy(byvars = ptCellVars,dt = cells, suffix = 'all')
  
  mkProportionBy(
    byvars = c(ptCellVars, 'is_epithelial'),
    dt = cells, suffix = 'isEpi')
  
  proportionVars <- grep('proportion_', names(cells), value = T)
  phenotypeVars <- grep('meta_|phenotype|colours', names(cells), value = T)
  indicatorVars <- grep('^is_', names(cells), value = T)
  indicatorVars <- setdiff(indicatorVars, c('is_normal', 'is_dcis'))
  
  toKeep <- c(ptCellVars, proportionVars, phenotypeVars, indicatorVars)
  proportionsOut <- cells[, .SD, .SDcols = toKeep]
  proportionsOut <- melt(proportionsOut,
                         id.vars = c('ImageNumber', phenotypeVars, indicatorVars),
                         measure.vars = proportionVars,
                         value = 'proportion',
                         variable = 'type'
  )	
  proportionsOut[, type := gsub('proportion_', '', type)]
  proportionsOut[type == 'isEpi', 
                 type := ifelse(is_epithelial, 'tumour', 'stroma')]
  proportionsOut <- proportionsOut[, .SD[1], by = .(ImageNumber, meta_description, type)]
  proportionsOut[, type := gsub(' ', '', type)]
  
  proportionsOut[, check := sum(proportion), by = .(ImageNumber, type)]
  stopifnot(all.equal(rep(1, nrow(proportionsOut)), proportionsOut[['check']]))
  proportionsOut[, check := NULL]
  
  outfile <- here('scratch/cellPhenotypeProportions.fst')
  write_fst(proportionsOut, outfile)
  
  # Cell phenotype proportions
  cellPhenotypes <- here('scratch/cellPhenotypeProportions.fst')
  cellPhenotypes <- read_fst(cellPhenotypes, as.data.table = T)
  cellPhenotypes <- cellPhenotypes[grep('tumour|stroma', type)]
  cellPhenotypes <- cellPhenotypes[, .(ImageNumber, type, meta_description, proportion)]
  allCombinations <- adt(expand.grid(unique(cellPhenotypes[,ImageNumber]), unique(cellPhenotypes[,meta_description])))
  setnames(allCombinations, c('ImageNumber', 'meta_description'))
  allCombinations <- merge(x = allCombinations, 
                           y = cellPhenotypes[!duplicated(meta_description), .(meta_description, type)], 
                           by = 'meta_description') 
  cellPhenotypes <- merge(x = allCombinations, y = cellPhenotypes, 
                          by = c('ImageNumber', 'meta_description', 'type'), all.x = T) 
  cellPhenotypes[, meta_description := gsub("[^A-Za-z0-9+\\-]", "", meta_description)]
  cellPhenotypes[, meta_description := paste0(meta_description, '_originalType')]
  cPh_tumour <- cellPhenotypes[type == 'tumour', unique(meta_description)]
  cPh_tme <- cellPhenotypes[type == 'stroma', unique(meta_description)]
  cellPhenotypes[is.na(proportion), proportion := 0]
  cellPhenotypes[, checkBothCompartments := sum(proportion), by = .(ImageNumber, type)]
  cellPhenotypes <- cellPhenotypes[checkBothCompartments > 0][, checkBothCompartments := NULL]
  cellPhenotypes <- dcast(cellPhenotypes, ImageNumber ~ meta_description, value.var = 'proportion')
  cellPhenotypes <- na.omit(cellPhenotypes) # samples that contain both tumour and stromal cells
  predictors <- merge(x = predictors, y = cellPhenotypes, by = 'ImageNumber')
  saveRDS(predictors, file = here('scratch/predictors_intfile.rds'))
  
  
  
  
  
  return(predictors)
  
}


apply_filtering <- function(predictors, percentage){

  generate_matrix<- function(df, subselection=NULL, NA_percentage=0){
    m <- as.data.frame(df)
    rownames(m) <- m$tnumber
    m <- m[,-c(1)]
    # filtering: with percentages`
    m <- m[, which(colMeans(!is.na(m)) > NA_percentage)]
    # Impute NANs before or after scaling?
    # m <- m %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
    
    if (!is.null(subselection)){
      m <- m %>% dplyr::select(all_of(subselection))
    }
    
    m <- scale(m)
    m <- replace(m,is.na(m),0)
    
    m <- as.data.frame(m)
    
    return(m)
  }
  
  shape_features <- getShapeFeaturesAlternative()
  scale_features <- getScaleFeaturesAlternative()
  estimations <- read_rds(here('scratch/estimationCountsAlternative.rds'))
  select_features <- rownames(estimations %>% filter(rank_percentage <= percentage))
  
  
  shape_features <- generate_matrix(shape_features,select_features, 0)
  scale_features <- generate_matrix(scale_features,select_features, 0)
  print(ncol(shape_features))
  
  shape_features$ImageNumber <- rownames(shape_features)
  scale_features$ImageNumber <- rownames(scale_features)
  
  shapePredictors <- setdiff(colnames(shape_features), c('metabric_id', 'ImageNumber'))
  scalePredictors <- setdiff(colnames(scale_features), c('metabric_id', 'ImageNumber'))
  
  shape_features <- shape_features %>%
    rename_with( .fn = function(.x){paste0("shape_", .x)},
                 .cols=all_of(shapePredictors))
  
  scale_features <- scale_features %>%
    rename_with( .fn = function(.x){paste0("scale_", .x)},
                 .cols=all_of(scalePredictors))
  
  shapePredictors <- setdiff(colnames(shape_features), c('metabric_id', 'ImageNumber'))
  scalePredictors <- setdiff(colnames(scale_features), c('metabric_id', 'ImageNumber'))
  
  shape_and_scale <- merge(shape_features, scale_features, by=c('ImageNumber'))
  shape_and_scale$ImageNumber = as.numeric(shape_and_scale$ImageNumber)
  
  predictors <- merge(x = predictors, y = shape_and_scale, by = 'ImageNumber')
  
  # Response intclusts
  intClust <- as.data.table(clinical_data)[, .(ImageNumber, IntClust)]
  intClust <- intClust[!is.na(IntClust)]
  for(i in sortunique(intClust$IntClust)) intClust[, eval(i) := as.numeric(IntClust == i)]
  intClust[, IntClust := NULL]
  
  # Response ER/Her2
  ERHer2 <- as.data.table(clinical_data)[, .(ImageNumber, ER_HER2_status)]
  ERHer2 <- ERHer2[!is.na(ER_HER2_status)]
  for(i in sortunique(ERHer2$ER_HER2_status)) ERHer2[, eval(i) := as.numeric(ER_HER2_status == i)]
  ERHer2[, ER_HER2_status := NULL]
  
  
  # Response PAM groups
  PAM50 <- as.data.table(clinical_data)[, .(ImageNumber, PAM50)]
  PAM50 <- PAM50[!is.na(PAM50)]
  for(i in sortunique(PAM50$PAM50)) PAM50[, eval(i) := as.numeric(PAM50 == i)]
  PAM50[, PAM50 := NULL]
  
  toModel <- merge(x = predictors, y = ERHer2, by = 'ImageNumber', all.x = T)
  toModel <- merge(x = toModel, y = PAM50, by = 'ImageNumber', all.x = T)
  saveRDS(predictors, file=here('scratch/model_input.rds'))
  
  return(toModel)
}


getPredictorName <- function(input_predictors){
  predictor_name <- 'unknown'
  
  if (length(setdiff(input_predictors, cPh_tme))==0){
    predictor_name <- 'TME_proportion'
  }
  if (length(setdiff(input_predictors, cPh_tumour))==0){
    predictor_name <- 'tumour_proportion'
  }
  if (length(setdiff(input_predictors, cPh_tme_original))==0){
    predictor_name <- 'originalTME_proportion'
  }
  if (length(setdiff(input_predictors, cPh_tumour_original))==0){
    predictor_name <- 'originalTumour_proportion'
  }
  
  if (length(setdiff(input_predictors, TMEStructurePredictors))==0){
    predictor_name <- 'TME_structures'
  }
  if (length(setdiff(input_predictors, nwPredictorsTME))==0){
    predictor_name <- 'TME_nw'
  }
  if (length(setdiff(input_predictors, nwPredictorsTumour))==0){
    predictor_name <- 'tumour_nw'
  }
  if (length(setdiff(input_predictors, shapePredictors))==0){
    predictor_name <- 'shape'
  }
  if (length(setdiff(input_predictors, scalePredictors))==0){
    predictor_name <- 'scale'
  }
  if (length(setdiff(c(shapePredictors, scalePredictors), input_predictors))==0){
    predictor_name <- 'shape_and_scale'
  }
  if (length(setdiff(input_predictors, cPh_tme))>0 & length(setdiff(input_predictors, cPh_tme)) < length(input_predictors)){
    predictor_name <- 'combination'
  }

  return(predictor_name)
}

fitModel <- function(response, predictors, Weights, dt){
  dt <- dt[, .SD, .SDcols = c(response, predictors, Weights)]	
  dt <- na.omit(dt)
  weights <- dt[, get(Weights)]
  response <- as.numeric(dt[, get(response)])
  predictors <- dt[, .SD, .SDcols = predictors]
  lassoFit <- cva.glmnet(as.matrix(predictors), as.matrix(response), weights = weights, 
                         family = 'binomial', nlambda = 100, nfolds = 10)
  return(lassoFit)
}

getAUC <- function(
    responseVar, # character column name
    predictors, # character vector of column names
    Weights, # character column name
    testTrain # list of dts named 'train' and 'test' 
){
  
  sampleSizeTrain <- testTrain$train[!is.na(get(responseVar)), .N]
  samplePositiveTrain <- testTrain$train[get(responseVar) == 1, .N]
  
  sampleSizeTest <- testTrain$test[!is.na(get(responseVar)), .N]
  samplePositiveTest <- testTrain$test[get(responseVar) == 1, .N]
  
  modelFit <- fitModel(response = responseVar, predictors = predictors, 
                       Weights = Weights, dt = testTrain$train)
  
  
  
  # Select alpha with best performance
  min_cvm <- min(modelFit[["modlist"]][[1]][["cvm"]])
  a <- modelFit$alpha[1]
  index <- 1
  for (m in 2:length(modelFit$alpha)){
    pot_min_cvm <- min(modelFit[["modlist"]][[m]][["cvm"]])
    if (pot_min_cvm < min_cvm){
      min_cvm <- pot_min_cvm
      a <- modelFit$alpha[[m]]
      index <- m
    }
  }

  predictor_name <- getPredictorName(predictors)
  saveRDS(predictor_name, here('scratch/predictorName.rds'))

  # coefList <- coef(modelFit, s=coefList$lambda.1se) #For cv.glmnet
  coefList <- coef(modelFit[['modlist']][[index]], s = "lambda.min") #For cva.glmnet
  coefList <- data.frame(coefList@Dimnames[[1]][coefList@i+1],coefList@x)
  names(coefList) <- c('var','val')
  coefList <- coefList %>% mutate(response = responseVar)
  coefList <- coefList %>% mutate(predictor = predictor_name)
  
  

	# Save files during final prediction
  save_rds_archive(coefList, file=paste(here('scratch/coeflists_run2/coefList_'), gsub('+', 'plus', responseVar, fixed=T), '_', predictor_name, '.rds', sep=''),last_modified = T, with_time = T)
  
  testDt <- testTrain$test[, .SD, .SDcols = c(responseVar, predictors)]
  
  testDt <- na.omit(testDt)
  testMat <- as.matrix(testDt[, .SD, .SDcols = predictors])
  predictions <- predict(modelFit, testMat, type = 'response', s = 'lambda.min',alpha=a) #For cva.glmnet
  # predictions <- predict(modelFit, testMat, type = 'response', s = 'lambda.min') #For cv.glmnet
  trueLabels <- as.numeric(testDt[, get(responseVar)])
  predictObj <- prediction(predictions, trueLabels)
  auc <- performance(predictObj,"auc") 
  auc <- as.numeric(auc@y.values)
  
  # roc <- pROC::roc(response = trueLabels, predictor = predictions, direction= "<")
  
  # Save files during filtering iterations
  # save_rds_archive(roc, file= paste(here('scratch/ROCS_run2/'),predictor_name ,'/roc_filter',length(predictors), '_' ,gsub('+', 'plus', responseVar, fixed=T), '_', predictor_name, '.rds', sep=''), last_modified = T, with_time = T)
  
  # Save files during test/train cross-validation
  # save_rds_archive(roc, file= paste(here('scratch/ROCS_run2/roc_'),responseVar, '_', predictor_name, '.rds', sep=''), last_modified = T, with_time = T)
  
  # Save files during final prediction
  # saveRDS(roc, file= paste(here('scratch/ROCS/roc_'),responseVar, '_', predictor_name, '.rds', sep=''))
  
  return(data.table(
    sampleSizeTrain = sampleSizeTrain, samplePositiveTrain = samplePositiveTrain,
    sampleSizeTest = sampleSizeTest, samplePositiveTest = samplePositiveTest,
    response = responseVar, auc = auc, alpha = a))
}


mkAUCsTable <- function(Rep, predictors, toModel){
  mkAUC <- function(responseVar, predictors, Weights, testTrain){
    tryCatch(
      expr = {getAUC(responseVar, predictors, Weights, testTrain)},
      error = function(e){return(e)}
    )
  }
  
  AUCs <- rbindlist(mclapply(responseVars, mkAUC, predictors = predictors, Weights = 'weights', testTrain = toModel, mc.cores = 12))
  AUCs[, rep := Rep]
  return(AUCs)
}


# Fit modeltypes
fitAndEstimate <- function(selection_predictors, name_predictors, outfile, df){
  RNGkind("L'Ecuyer-CMRG")
  mc.reset.stream()
  doAUCs <- function(predictors){
    Reps <- seq_len(1)
    AUCs <- rbindlist(lapply(Reps, mkAUCsTable, predictors = predictors, toModel=df))
    return(AUCs)
  }
  
  AUCs<- lapply(selection_predictors, doAUCs)
  names(AUCs) <- name_predictors
  
  modeltypes <- seq_len(length(AUCs))
  names(modeltypes) <- names(AUCs)
  AUCs <- rbindlist(lapply(names(AUCs), function(name){AUCs[[name]][, type := name]}))
  
  AUCs[, facet_by := modeltypes[type]]
  AUCs[, facet_by := factor(facet_by, 
                            levels = modeltypes, labels = names(modeltypes), ordered = TRUE)]
  write_fst(AUCs, outfile)
  
}

### Visualization methods ###


collect_intclust <- function(AUCs){
  AUCs <- read_fst(AUCs)
  
  intClustAUCs <- AUCs[grep('IntClust', response)]
  # Plot ic10 AUCs
  intClustCols <- getSubtypeCols()
  intClustCols <- intClustCols$IntClust
  intClustAUCs[, yaxis := as.numeric(gsub('IntClust ', '', response))]
  intClustAUCs[yaxis > 4, yaxis := yaxis + 2]
  intClustAUCs[grep('4\\+', response), yaxis := 4]
  intClustAUCs[grep('4\\-', response), yaxis := 5]
  intClustAUCs[grep('5\\+', response), yaxis := 6]
  intClustAUCs[grep('5\\-', response), yaxis := 7]
  intClustAUCs <- intClustAUCs[order(type, yaxis),]
  
  intClustAUCs[, yaxis := factor(yaxis, 
                                 levels = unique(yaxis), labels = unique(response), ordered = TRUE)]
  intClustAUCs[, colours := intClustCols[response]]
  intClustAUCs[, yaxis := reorder(yaxis, dplyr::desc(yaxis))]
  
  x <- c("Cell phenotype (Tumour)","Cell phenotype (TME)","TME Structures", "Network properties (Tumour)", "Network properties (TME)","shape features","scale features", "shape and scale features",  "Combined")
  
  intClustAUCs <- intClustAUCs %>%
    arrange(sapply(type, function(y) which(y == x)))
  
  return(intClustAUCs)
  
}

collect_pam <- function(AUCs){
  AUCs <- read_fst(AUCs)
  
  intClustAUCs <- AUCs %>% filter(response %in% c('Luminal A','Luminal B', 'HER2', 'Basal', 'Normal-like'))
  intClustCols <- getSubtypeCols()
  intClustCols <- intClustCols$PAM50
  
  intClustAUCs <- intClustAUCs %>% mutate(yaxis = ifelse(response=='Basal', 1, ifelse(response=='HER2', 2,ifelse(response=='Luminal A', 3, ifelse(response=='Luminal B', 4, 5)))))
  
  intClustAUCs <- intClustAUCs %>% mutate(yaxis = factor(yaxis, 
                                                         levels = unique(yaxis), labels = unique(response), ordered = TRUE))
  intClustAUCs <- intClustAUCs %>% mutate(colours = intClustCols[response])
  intClustAUCs <- intClustAUCs %>% dplyr::arrange(desc(yaxis))
  
  x <- c("Cell phenotype (Tumour)","Cell phenotype (TME)","TME Structures", "Network properties (Tumour)", "Network properties (TME)","shape features","scale features", "shape and scale features",  "Combined")
  
  intClustAUCs <- intClustAUCs %>%
    arrange(sapply(type, function(y) which(y == x)))
  
  return(intClustAUCs)
}

collect_MolecularSubtypes <- function(AUCs){
  AUCs <- read_fst(AUCs)
  
  intClustAUCs <- AUCs %>% filter(response %in% c('ER+HER2+','ER-HER2+','ER+HER2-','ER-HER2-'))
  intClustCols <- getSubtypeCols()
  intClustCols <- intClustCols$molecularSubtypes
  
  intClustAUCs <- intClustAUCs %>% mutate(yaxis = ifelse(response=='ER+HER2+', 1, ifelse(response=='ER-HER2+', 2,ifelse(response=='ER+HER2-', 3, 4))))
  
  intClustAUCs <- intClustAUCs %>% mutate(yaxis = factor(yaxis, 
                                                         levels = unique(yaxis), labels = unique(response), ordered = TRUE))
  intClustAUCs <- intClustAUCs %>% mutate(colours = intClustCols[response])
  intClustAUCs <- intClustAUCs %>% dplyr::arrange(desc(yaxis))
  
  x <- c("Cell phenotype (Tumour)","Cell phenotype (TME)","TME Structures", "Network properties (Tumour)", "Network properties (TME)","shape features","scale features", "shape and scale features",  "Combined")
  
  intClustAUCs <- intClustAUCs %>%
    arrange(sapply(type, function(y) which(y == x)))
  
  return((intClustAUCs))
  
}


generate_plot <- function(df){
  aucPlot <- ggplot() +
    geom_vline(xintercept = seq(0.4, 0.8, 0.1), colour = 'lightgrey', size = 0.3) +
    geom_vline(xintercept = 0.5, colour = 'steelblue', size = 0.5, linetype = 'dotted') +
    geom_point(data = df, aes(x = auc, y = yaxis, fill=colours),alpha=1, size = 3, pch = 21) +
    xlim(0.2,1) +
    scale_fill_identity() +
    scale_colour_identity() +
    theme_prj(base_line_size = 0.25) +
    theme(legend.position='left',
          panel.background = element_rect(colour = 'black', size = 0.25),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.x = unit(1,'mm'),
          plot.margin = unit(c(0,0.25,0,0), 'mm')) +
    facet_wrap(vars(facet_by), nrow = 1) +
    labs(x = 'AUC')
  
  
  mkCountBar <- function(DaT){
    barCount <- ggplot(data = DaT) +
      geom_bar( 
        aes(y = samplePositiveTest, x = yaxis), 
        width = 0.8, stat = 'identity', fill = 'grey') +
      geom_text(aes(y = samplePositiveTest, x = yaxis, label = samplePositiveTest), 
                hjust = 0, nudge_y = 1) +
      facet_wrap(vars(facet_by)) +
      theme_prj() +
      theme(strip.text = element_blank(),
            plot.margin = unit(c(0,5,0,0), 'mm'),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),) +
      labs(y = bquote(italic('n')~'tumours')) +
      coord_flip(clip = 'off')
    return(barCount)	
  }
  
  df <- df[order(df$response), ]
  df <- df[!duplicated(df$response), ]
  
  
  barCount <- mkCountBar(df)
  
  icOut <- plot_grid(aucPlot, barCount, align = 'h', rel_widths = c(9, 0.4))
  
  return(icOut)
  
}

generate_plot_alphas <- function(df){
  aucPlot <- ggplot() +
    geom_vline(xintercept = seq(0.4, 0.8, 0.1), colour = 'lightgrey', size = 0.3) +
    geom_vline(xintercept = 0.5, colour = 'steelblue', size = 0.5, linetype = 'dotted') +
    geom_point(data = df, 
               aes(x = auc, y = yaxis, fill=alpha),alpha=0.5, size = 2.5, pch = 21) +
    scale_fill_gradient2(midpoint=0.5, low="blue", mid="white",
                         high="red", space ="Lab", name='alpha' ) + 
    theme_prj(base_line_size = 0.25) +
    theme(legend.position='left',
          panel.background = element_rect(colour = 'black', size = 0.25),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.x = unit(1,'mm'),
          plot.margin = unit(c(0,0.25,0,0), 'mm')) +
    facet_wrap(vars(facet_by), nrow = 1) +
    labs(x = 'AUC')
  
  
  mkCountBar <- function(DaT){
    barCount <- ggplot(data = DaT) +
      geom_bar( 
        aes(y = samplePositiveTest, x = yaxis), 
        width = 0.8, stat = 'identity', fill = 'grey') +
      geom_text(aes(y = samplePositiveTest, x = yaxis, label = samplePositiveTest), 
                hjust = 0, nudge_y = 1) +
      facet_wrap(vars(facet_by)) +
      theme_prj() +
      theme(strip.text = element_blank(),
            plot.margin = unit(c(0,5,0,0), 'mm'),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),) +
      labs(y = bquote(italic('n')~'tumours')) +
      coord_flip(clip = 'off')
    return(barCount)	
  }
  df <- df[order(df$response), ]
  df <- df[!duplicated(df$response), ]
  
  
  barCount <- mkCountBar(df)
  
  icOut <- plot_grid(aucPlot, barCount, align = 'h', rel_widths = c(9, 0.4))
  
  return(icOut)
  
}

generate_plot_filtering <- function(df){
  aucPlot <- ggplot() +
    geom_vline(xintercept = seq(0.4, 0.8, 0.1), colour = 'lightgrey', size = 0.3) +
    geom_vline(xintercept = 0.5, colour = 'steelblue', size = 0.5, linetype = 'dotted') +
    geom_point(data = df, 
               aes(x = auc, y = yaxis, fill=colours), size = 2.5, pch=21) +
    scale_fill_gradient2(low="white",
                         high="black", space ="Lab", name='% features included' ) + 
    theme_prj(base_line_size = 0.25) +
    theme(legend.position='left',
          panel.background = element_rect(colour = 'black', size = 0.25),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.x = unit(1,'mm'),
          plot.margin = unit(c(0,0.25,0,0), 'mm')) +
    facet_wrap(vars(facet_by), nrow = 1) +
    labs(x = 'AUC')
  
  
  mkCountBar <- function(DaT){
    barCount <- ggplot(data = DaT) +
      geom_bar( 
        aes(y = samplePositiveTest, x = yaxis), 
        width = 0.8, stat = 'identity', fill = 'grey') +
      geom_text(aes(y = samplePositiveTest, x = yaxis, label = samplePositiveTest), 
                hjust = 0, nudge_y = 1) +
      facet_wrap(vars(facet_by)) +
      theme_prj() +
      theme(strip.text = element_blank(),
            plot.margin = unit(c(0,5,0,0), 'mm'),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),) +
      labs(y = bquote(italic('n')~'tumours')) +
      coord_flip(clip = 'off')
    return(barCount)	
  }
  
  df <- df[order(df$response), ]
  df <- df[!duplicated(df$response), ]
  
  barCount <- mkCountBar(df)
  
  icOut <- plot_grid(aucPlot, barCount, align = 'h', rel_widths = c(3, 0.4))
  
  return(icOut)
  
}



generate_plot_withboxplot <- function(df){
  confidenceIntervals <- df %>% group_by(yaxis, facet_by) %>% summarise(mean_auc = mean(auc), sd_auc = sd(auc)) %>% mutate(error = qnorm(0.975)*sd_auc/sqrt(5)) %>% mutate(up = mean_auc + error) %>% mutate(down = mean_auc - sd_auc)
  
  aucPlot <- ggplot() +
    geom_vline(xintercept = seq(0.4, 0.8, 0.1), colour = 'lightgrey', size = 0.3) +
    geom_vline(xintercept = 0.5, colour = 'steelblue', size = 0.5, linetype = 'dotted') +
    geom_point(data=confidenceIntervals, aes(x=up, y=yaxis), pch=3) +
    geom_point(data=confidenceIntervals, aes(x=down, y=yaxis), pch=3) +
    geom_segment(data=confidenceIntervals, aes(x=down,xend=up, y=yaxis, yend=yaxis))  + 
    geom_point(data = df %>% filter(type == 'test_set'), aes(x = auc, y = yaxis, fill=colours),alpha=1, size = 3, pch = 21) +
    geom_point(data = confidenceIntervals, aes(x = mean_auc, y = yaxis),colour='black',alpha=1, size = 3, pch = 3) +
    xlim(0.2,1) +
    theme_prj(base_line_size = 0.25) +
    theme(legend.position='none',
          panel.background = element_rect(colour = 'black', size = 0.25),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.x = unit(1,'mm'),
          plot.margin = unit(c(0,0.25,0,0), 'mm')) +
    facet_wrap(vars(facet_by), nrow = 1) +
    labs(x = 'AUC') 
  
  mkCountBar <- function(DaT){
    barCount <- ggplot(data = DaT) +
      geom_bar( 
        aes(y = samplePositiveTest, x = yaxis), 
        width = 0.8, stat = 'identity', fill = 'grey') +
      geom_text(aes(y = samplePositiveTest, x = yaxis, label = samplePositiveTest), 
                hjust = 0, nudge_y = 1) +
      facet_wrap(vars(facet_by)) +
      theme_prj() +
      theme(strip.text = element_blank(),
            plot.margin = unit(c(0,5,0,0), 'mm'),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),) +
      labs(y = bquote(italic('n')~'tumours')) +
      coord_flip(clip = 'off')
    return(barCount)	
  }
  
  df <- df[order(df$response), ]
  df <- df[!duplicated(df$response), ]
  
  barCount <- mkCountBar(df)
  
  
  icOut <- plot_grid(aucPlot, barCount, align = 'h', rel_widths = c(9, 0.4))
  
  return(icOut)
  
}
