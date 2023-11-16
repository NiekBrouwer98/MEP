# This file contains the code to generate clean feature sets from the raw data output of the distance method
library(tidyverse)
source(here("MEP_UtilityFunctions.R"))


get_parameters <- function(path){
  params <- c()
  files <- list.files(path)
  for (f in files){
    tryCatch({
      success_model <- readRDS(paste(path,f,sep=''))
      params <- append(params, list(success_model[[3]]))
      name <- gsub('success_models_','',f)
      name <- gsub('.rds','',name)
      names(params)[length(params)] <- name
    },error=function(e){print(e)}, warning=function(w){print(w)})
    
  }
  params_result <- bind_rows(params, .id='phenotype_combo')
  return(params_result)
}


#' The first parameter estimation was done with thresholds [0,20,50,70,100]
#' The second parameter estimation was done with thresholds [0,2,5,10,20]
#' Merge the two estimations
#' @return Save parameter file in scratch dir
CompileALLParameters <- function(save=TRUE){

  all_parameters_firstrun <- get_parameters(here('scratch/success_models/')) %>%
    separate(phenotype_combo, into=c('phenotype_from','phenotype_to'),sep='_to_', remove = FALSE) %>%
    filter(a > 0.5) %>%
    filter(b > 8) %>%
    rename(shape =a , scale = b) %>%
    unite('unique_sample', c(tnumber, phenotype_combo), remove=F)
  
  all_parameters_secondrun <- get_parameters(here('scratch/success_models_secondrun/')) %>%
    separate(phenotype_combo, into=c('phenotype_from','phenotype_to'),sep='_to_', remove = FALSE) %>%
    filter(a > 0.5) %>%
    filter(b > 8) %>%
    rename(shape =a , scale = b) %>%
    unite('unique_sample', c(tnumber, phenotype_combo), remove=F)
  
  #Bind both results
  all_parameters <- rbind(all_parameters_secondrun, all_parameters_firstrun %>% filter(! unique_sample %in% all_parameters_secondrun$unique_sample))
  if(save){
    saveRDS(all_parameters, file = here('scratch/all_parameters.rds'))
  }
  
  return(all_parameters)

}


CompileALLParameters_run2 <- function(save=TRUE){
  all_parameters <- get_parameters(here('scratch/success_models_run2/')) %>%
    separate(phenotype_combo, into=c('phenotype_from','phenotype_to'),sep='_to_', remove = FALSE) %>%
    filter(a > 0.5) %>%
    filter(b > 8) %>%
    dplyr::rename(shape =a , scale = b) %>%
    unite('unique_sample', c(tnumber, phenotype_combo), remove=F)
  
  if(save){
    saveRDS(all_parameters, file = here('scratch/all_parameters_run2.rds'))
  }
  
  return(all_parameters)
  
}

#' Generate matrix from the input parameters
#' Keep in mind whether to include outliers
#' @return Save parameter file in scratch dir
transform_parameters_to_matrix <- function(all_parameters, shapeOutput_path, scaleOutput_path){
  shape_parameters <- tibble(tnumber = unique(all_parameters$tnumber))
  scale_parameters <- tibble(tnumber = unique(all_parameters$tnumber))
  
  for (c in unique(all_parameters$phenotype_combo)){
    shape_parameters <- left_join(x = shape_parameters, y= all_parameters %>% 
                                    dplyr::select(c(tnumber, phenotype_combo, shape)) %>% 
                                    filter(phenotype_combo == c), by='tnumber') %>% dplyr::select(-c(phenotype_combo))
    names(shape_parameters)[names(shape_parameters) == 'shape'] = c
    
    scale_parameters <- left_join(x = scale_parameters, y= all_parameters %>% 
                                    dplyr::select(c(tnumber, phenotype_combo, scale)) %>% 
                                    filter(phenotype_combo == c), by='tnumber') %>% dplyr::select(-c(phenotype_combo))
    names(scale_parameters)[names(scale_parameters) == 'scale'] = c
  }
  
  saveRDS(shape_parameters, shapeOutput_path)
  saveRDS(scale_parameters, scaleOutput_path)
}

getInitialParameters <- function(){
  combinations <- getCombinations()
  initial_parameters <- tibble()
  for (c in combinations){
      init_params <- readRDS(here(paste('scratch/init_parameters/init_params_', c,'.rds', sep='')))
      initial_parameters <- bind_rows(initial_parameters, init_params)
  }

  initial_parameters <- reshape2::dcast(initial_parameters, tnumber + combo ~ term, value.var = 'estimate')

  initial_parameters <- initial_parameters %>%
                        unite('unique_sample', c(tnumber, combo), remove=FALSE) %>%
                        separate(combo, into=c('phenotype_from','phenotype_to'),sep='_to_', remove=FALSE) %>%
                        dplyr::rename(phenotype_combo = combo)

  saveRDS(initial_parameters, file=here('scratch/initial_parameters.rds'))
  
  
}

getInitialParametersAlternative <- function(){
  combinations <- getCombinationsAlternative()
  initial_parameters <- tibble()
  for (c in combinations){
    init_params <- readRDS(here(paste('scratch/init_parameters_run2/init_params_', c,'.rds', sep='')))
    initial_parameters <- bind_rows(initial_parameters, init_params)
  }
  
  initial_parameters <- reshape2::dcast(initial_parameters, tnumber + combo ~ term, value.var = 'estimate')
  
  initial_parameters <- initial_parameters %>%
    unite('unique_sample', c(tnumber, combo), remove=FALSE) %>%
    separate(combo, into=c('phenotype_from','phenotype_to'),sep='_to_', remove=FALSE) %>%
    dplyr::rename(phenotype_combo = combo)
  
  saveRDS(initial_parameters, file=here('scratch/initial_parameters_run2.rds'))
  
  
}

