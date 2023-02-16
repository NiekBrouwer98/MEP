# Required packages for code
library(tidyverse)
library(ggpubr)
library(readxl)
library(reshape2)
library(dbscan)
library(broom)
library(zoo)
library(knitr)
library(magrittr)
library(NMF,quietly = T)
library(MASS)
library(car)
library(stats)
library(foreach)
library(spatstat)
# library(raster)
library(nlme)
library(dplyr)
library(gridExtra)
library(here)
library(tictoc)

# Required packages for parallelization
library(doMC)

cores = detectCores()
registerDoMC(cores[1])

source(here("UtilityFunctions.R"))
library(fst)
library(data.table)
library(ggplot2)
library(assertthat)


initial_params <- readRDS(here('scratch/initial_parameters.rds'))
n_cells <- readRDS(here('scratch/n_cells.rds'))
xy_data_exc_oct <- readRDS(here('scratch/xy_data_exc_oct.rds'))

# Define phenotype combinations to study in data
# combinations <- unique(xy_data_exc_oct %>% pull(phenotype_combo))
combinations <- readRDS(here('scratch/combinations_selection.rds'))
combinations_selection <- combinations
# combinations_selection <- c("CK8-18^{+} ER^{hi}_to_CK8-18^{+} ER^{hi}")

threshold_pixels <- 300
threshold_microns <- threshold_pixels

# Define maximum Shape (A_max) and maximum Scale (B_max) parameters for your fits
B_max <- 500 # max scale parameter (0,500]
A_max <- 10  # max shape parameter (0,10]

# Define Weibull and parameterized Weibull functions, and functions to convert between parameterized and
## non parameterized Weibull
## (Parameterization was done to ease the fitting procedure, in a non-restricted manner)
convert_params <- function(a,b){
  A=-log(A_max/a - 1)
  A = as.numeric(A)
  B=-log(B_max/b-1)
  B = as.numeric(B)
  return(c(A=A,B=B))
}

convert_params_reverse <- function(A,B){
  a=A_max/(1+exp(-A))
  b=B_max/(1+exp(-B))
  return(c(a=a,b=b))
}
# Define Weibull distribution function in a parameterized way:
nform2 <- ~(( ( A_max/(1+exp(-A)) )  / ( B_max/(1+exp(-B)) )) * (( x/( B_max/(1+exp(-B))) )^(( A_max/(1+exp(-A)) )-1)) * exp(-1* (x/( B_max/(1+exp(-B))))^( A_max/(1+exp(-A)) )))

weibfunc2 <- function(x,A,B){
  a = A_max/(1+exp(-A))
  b = B_max/(1+exp(-B)) 
  return(((a/b) * ((x/b)^(a-1)) * exp(- (x/b)^a)))
}

nfun2 <- deriv(nform2,namevec=c("A","B"),
               function.arg=c("x","A","B"))


# Define function to fit a Weibull model to x/y coordinates of normalized 1-NN distance
## distribution coordinates, filtering out samples with n_cells < threshold_filter_cells
#' @param xy_coordinates_curves 1-NN distribution coordinate dataframe
#' @param combi Cell combination of interest: 'CellFROM_to_CellTO'
#' @param n_cells_df Dataframe with number of cells per sample (for each cell type)
#' @param threshold_filter_cells Threshold to filter out samples depending on the cell number
#' @param initial_params_df Dataframe with initial Weibull parameters
fit_data <- function(xy_coordinates_curves, combi, n_cells_df, threshold_filter_cells,
                     initial_params_df){
  message(combi)
  xy_data_filt <- xy_coordinates_curves %>% filter(phenotype_combo == combi) %>% as_tibble
  
  xy_data_filt <- xy_data_filt %>% 
    left_join(n_cells_df %>% mutate(n_from=n) %>% dplyr::select(-n), by=c('tnumber','phenotype_from'='phenotype')) %>%
    left_join(n_cells_df %>% mutate(n_to=n) %>% dplyr::select(-n), by=c('tnumber','phenotype_to'='phenotype'))
  
  threshold_cells <- threshold_filter_cells
  message(xy_data_filt %>% filter(n_from<=threshold_cells & n_to <= threshold_cells) %>%
            mutate(tnumber=paste(tnumber, 'nFROM', as.character(n_from), 'nTO',as.character(n_to))) %>% 
            pull(tnumber) %>% unique %>% paste(collapse=' / '))
  
  xy_data_filt <- xy_data_filt %>% filter(n_from>threshold_cells & n_to > threshold_cells)
  
  # Add missing observations in data and set to 0
  ## i.e. if no distances were observed beyond 100 microns, the normalized scaled counts 
  ## from 100 to 300 microns will be set to 0 
  # include 0,0 or only 1,0??
  xy_data_filt <- expand.grid(x=c(1,xy_data_filt %>%pull(distance_window) %>% unique), 
                              tnumber=xy_data_filt %>%pull(tnumber) %>% unique,
                              phenotype_combo=combi) %>%
    left_join(xy_data_filt, by=c('x','tnumber', 'phenotype_combo')) %>%
    mutate(y=ifelse(is.na(y),0,y))
  
  mean_params <- initial_params_df %>% filter(combo == combi) %>%
    filter(tnumber %in% unique(xy_data_filt$tnumber)) %>% 
    group_by(term) %>% dplyr::summarise(mean_estimate = mean(estimate))
  
  start_params <- c(a=mean_params %>% filter(term == 'shape') %>% pull(mean_estimate),
                    b=mean_params %>% filter(term == 'scale') %>% pull(mean_estimate))
  
  fit52_oct <- tryCatch({
    nlme(y~nfun2(x,A,B),
         fixed=A+B~1,
         random=list(tnumber=pdDiag(A+B~1)), # or pdDiag, pdLogChol
         data=xy_data_filt %>% mutate(tnumber = factor(tnumber)),
         start=convert_params(start_params['a'], start_params['b']),
         method='REML', control = nlmeControl(maxIter = 1000, 
                                              returnObject = F,
                                              msMaxIter=200,
                                              pnlsTol=1e-1, 
                                              tolerance=1e-6, opt="nlm"))
  },
  error=function(cond){
    message(paste('failed',combi))
    return(combi)
  })
  # }
  return(list(model=fit52_oct, xy_data=xy_data_filt))
  
}

tic('Weibull parameter estimation')

success_models <- list()

# Dataframe to store fitted Weibull coefficients
weib_coefs <- data.frame(matrix(nrow=0,ncol=4))
colnames(weib_coefs) <- c('a','b','sample','phenotype_combo')
weib_coefs$a %<>% as.numeric()  # shape
weib_coefs$b %<>% as.numeric()  # scale
weib_coefs$sample %<>% as.character()
weib_coefs$phenotype_combo %<>% as.character()

# Function to store x/y coordinates used for modelling
xy_data_weib <- xy_data_exc_oct %>% slice(0)
xy_data_weib <- xy_data_weib %>% add_column(yhat=as.numeric())


dir.create(here('scratch/success_models'))

success_models <- mclapply(combinations_selection, mc.cores = detectCores(), function(combi){
  for(threshold_cells in c(0,20,50,70,100)){
      message(threshold_cells)
      pipeline_run <- fit_data(xy_data_exc_oct, combi, n_cells,
                               threshold_cells, initial_params)

      fitted_model <- pipeline_run[['model']]

      if(!is(fitted_model, 'character')){
        new_data_weib_param <- as_tibble(ranef(fitted_model)) %>%
          mutate(sample=rownames(ranef(fitted_model))) %>%
          mutate(A=A + fixef(fitted_model)['A'],
                 B=B+ fixef(fitted_model)['B']) %>%
          mutate(phenotype_combo=combi) %>%
          mutate(a=10/(1+exp(-A)),b=500/(1+exp(-B)))


        # x/y coordinates of input / output model:
        xy_data_weib <- pipeline_run[['xy_data']] %>% mutate(yhat=fitted(fitted_model))
        
        # Process parameters
        fixefect <- fixef(pipeline_run[['model']])
        ranef <- ranef(pipeline_run[['model']])
        
        final_parameters <- ranef %>%
          mutate(tnumber= rownames(ranef),
                 phenotype_combo=combi) %>% 
          mutate(A = A + fixefect[['A']],
                 B = B + fixefect[['B']]) %>%
          mutate(a=10/(1+exp(-A)),b=500/(1+exp(-B)))
        
        result <- list(new_data_weib_param, xy_data_weib, final_parameters)
        # result <- list(result)
        # result <- setNames(result, combi)


        saveRDS(result, file = here(paste('scratch/success_models/success_models_', combi, '.rds', sep = '')))


        return(result)

      }
  }
  
  failed_combi <- list(list())
  # failed_combi <- list(failed_combi)
  # failed_combi <- setNames(failed_combi, combi)
  
  return(failed_combi)
})
    

names(success_models) <- combinations_selection
    
toc()

saveRDS(success_models, file = here('scratch/success_models.rds'))

#If we are here, all combinations are estimated and we don't have to save intermediate results
unlink(here("scratch/success_models"),recursive=TRUE)

