# Required packages for code
library(ComplexHeatmap)
library(tidyverse)
library(plyr)
library(cold) #function fixeff()
library(lme4)
library(glmmTMB)
library(fitdistrplus)
library(tidyr)
library(cold)
library(broom)
library(ggpubr) # stat_compare_means
library(here)
library(fst)
library(tictoc)

# Required packages for parallelization
library(doMC)

cores = detectCores()
registerDoMC(cores[1])

all_distances_data <- readRDS(file  = here('scratch/AUCScaledSlides_300_ALL.rds'))

all_distances_data <- all_distances_data %>% as_tibble %>%
  mutate(distance_window = WinMean) %>%
  mutate(phenotype_combo = paste(phenotype_from, phenotype_to, sep='_to_')) %>%
  dplyr::select(tnumber, phenotype_combo, `N.per.mm2.scaled`, distance_window) %>%
  mutate(new = `N.per.mm2.scaled` * 1000) %>%
  mutate(new =round(new))

combinations <- unique(all_distances_data %>% pull(phenotype_combo))
combinations_selection <- combinations
combinations_selection <- c("CK8-18^{+} ER^{hi}_to_CK8-18^{+} ER^{hi}")
# combinations_selection <- combinations[! combinations %in% c("Fibroblasts_to_CK^{+} CXCL12^{+}", "CK8-18^{hi}CXCL12^{hi}_to_CK^{+} CXCL12^{+}")]
saveRDS(combinations_selection, here('scratch/combinations_selection.rds')) #For same use in estimation

message(combinations_selection)

tic('parameter initialization')

all_combos_dists <- do.call(rbind, mclapply(combinations_selection, mc.cores = detectCores(), function(combi){
  # Recreate 1-NN histogram from coordinates (this is required for function "fitdistrplus")
  dists_onecombi <- tibble()
  message(combi)
  for(tn in unique(all_distances_data %>% pull(tnumber))){
    all_dists_tnum <- list()
    dists <- c(1,all_distances_data %>% filter(tnumber == tn) %>% filter(phenotype_combo == combi) %>% pull(distance_window))
    times <- c(1, all_distances_data %>% filter(tnumber == tn) %>% filter(phenotype_combo == combi) %>% pull(new))

    if(length(times) == 0){
      dists <- c(1:298)
      times <- rep(0,length(c(1:298)))
    }
    for(i in 1:length(times)){
      all_dists_tnum <- c(all_dists_tnum, rep(dists[i], times[i]))
    }
    aldists <- unlist(all_dists_tnum %>% as.numeric)
    aldistsdf <- as_tibble(data.frame(dists = aldists,
                                      tnumber = rep(tn, length(aldists)),
                                      combi = rep(combi, length(aldists))))
    dists_onecombi <- bind_rows(dists_onecombi, aldistsdf)
  }
  return(dists_onecombi)
  }))

saveRDS(all_combos_dists, file = here('scratch/all_combos_dists_one.rds'))
print('first step complete')

dir.create(here('scratch/init_parameters'))

initial_parameters <- do.call(rbind, mclapply(combinations_selection, mc.cores = detectCores(), function(c){
  message(c)
  initial_params <- data.frame(matrix(ncol=5, nrow=0))
  colnames(initial_params) <- c('term','estimate','std.error', 'tnumber','combo')
  initial_params <- initial_params %>% 
    mutate(term=as.character(term), estimate=as.numeric(estimate), `std.error`=as.numeric(`std.error`),tnumber=as.character(tnumber),
           combo=as.character(combo))
  
  for(tn in unique(all_combos_dists %>% pull(tnumber))){
    # Some vectors do not contain enough points to estimate a Weibull distribution
    # These vectors are not initialized
    tryCatch({
      params <- fitdist(all_combos_dists %>% filter(combi == c) %>% filter(tnumber == tn) %>% filter(dists < 100) %>% pull(dists), "weibull", lower = c(0, 0)) %>% summary
    class(params) <- c("fitdist", "fitdistr")
    params <- broom::tidy(params)
    params <- params %>% mutate(tnumber=rep(tn,nrow(params))) %>%
      mutate(combo=rep(c,nrow(params)))
    
    initial_params <- bind_rows(initial_params, params)
    
      },error=function(e){})
    
  }
  
  saveRDS(initial_params, here(paste('scratch/init_parameters/init_params_', c, '.rds', sep='' )))
  
  return(initial_params)
    
  }))

toc()
# Last run: 73876.134 s

saveRDS(initial_parameters, file = here('scratch/initial_parameters_one.rds'))

#If we are here, all combinations are initialized and we don't have to save intermediate results
unlink(here("scratch/init_parameters"),recursive=TRUE)

