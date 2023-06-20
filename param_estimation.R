# Required packages for code
library(tidyverse)
library(magrittr)
library(nlme)
library(dplyr)
library(here)
library(tictoc)

# Required packages for parallelization
library(doMC)

cores = detectCores()
registerDoMC(cores=8)

source(here("UtilityFunctions.R"))
source(here("MEP_UtilityFunctions.R"))

library(data.table)
library(fst)
library(R.utils)


###CONFIG
# structure_type <- 'Vascular stroma'
# structure_type <- 'Suppressed expansion' #REMOVE 'RUN2'
# structure_type <- 'APC enriched'
# structure_type <- 'Granulocyte enriched'
structure_type <- "TLSlike"

print(structure_type)

distance_data <- paste(here('scratch/AUCScaledSlides_300_'),gsub(' ', '_', structure_type, fixed=T) ,'.tsv', sep='')
initial_params_folder <- paste(here('scratch/init_parameters_'),gsub(' ', '_',structure_type,fixed=T), sep = '')
cells <- readRDS(paste(here('DATA/splittedCells_'), gsub(' ','_', structure_type,fixed = T), '.rds',sep='')) %>% dplyr::select(-c(ImageNumber)) %>% dplyr::rename(ImageNumber = Split_ImageNumber)
outfile <- paste(here('scratch/success_models_'),gsub(' ', '_',structure_type,fixed=T),sep='')

n_cells <- getCellCounts(cells) %>% dplyr::rename(phenotype = meta_description) %>% filter(n > 0)

slides <- read_tsv(distance_data)
##

# saveRDS(slides, file=here('scratch/AUCScaledSlides_300_ALL_run2.rds'))
xy_data <- slides %>%
  mutate(distance_window = WinMean) %>%
  mutate(x=distance_window, y=`N.per.mm2.scaled`)

# Define phenotype combinations to study in data
combinations <- readRDS(here('scratch/combinations_selection_run2.rds'))

# Directory for results
# Directory for intermediate results
if (!dir.exists(paste(here('scratch/success_models_'),gsub(' ', '_',structure_type,fixed=T), sep = ''))){
  dir.create(paste(here('scratch/success_models_'),gsub(' ', '_',structure_type,fixed=T), sep = ''))
}else{
  message("dir exists")
}

finished_files <- list.files(paste(here('scratch/success_models_'),gsub(' ', '_',structure_type,fixed=T),'/', sep = ''))
finished_combinations_int <-  finished_files %>% str_replace("success_models_", "")
finished_combinations_dir <- finished_combinations_int %>% str_replace(".rds", "")
finished_combinations <- unique(finished_combinations_dir)

combinations_selection <- combinations[combinations %in% finished_combinations == FALSE]

message(length(combinations_selection), " combinations: ", combinations_selection)

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
    res <- withTimeout({
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
    },timeout = 1800)
    },
  error=function(cond){
    message(paste('failed',combi))
    return(combi)
  }, TimeoutException = function(ex){
    message("Timemout")
    return(combi)
  })
  return(list(model=fit52_oct, xy_data=xy_data_filt))
  
}

tic('Weibull parameter estimation')

# Dataframe to store fitted Weibull coefficients
weib_coefs <- data.frame(matrix(nrow=0,ncol=4))
colnames(weib_coefs) <- c('a','b','sample','phenotype_combo')
weib_coefs$a %<>% as.numeric()  # shape
weib_coefs$b %<>% as.numeric()  # scale
weib_coefs$sample %<>% as.character()
weib_coefs$phenotype_combo %<>% as.character()

# Function to store x/y coordinates used for modelling
xy_data_weib <- xy_data %>% slice(0)
xy_data_weib <- xy_data_weib %>% add_column(yhat=as.numeric())

success_models <- mclapply(combinations_selection, mc.cores = 8, function(combi){
  # Skip combinations that take to long to avoid hold-up
  for(threshold_cells in c(0, 5,10,20,50)){
      message(threshold_cells)
      initial_params <- readRDS(paste(initial_params_folder, '/init_params_', combi, '.rds', sep=''))
      pipeline_run <- fit_data(xy_data, combi, n_cells,
                               threshold_cells, initial_params)

      fitted_model <- pipeline_run[['model']]

      # If fitted model succeeded: process
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

        #Save intermediate result
        saveRDS(result, file = paste(outfile, '/success_models_', combi, '.rds', sep = ''))

        return(list())
      }
      
  }
    
    # Parameter estimation failed
    saveRDS(list('FAILEDESTIMATION',list(),list()), file = paste(outfile, '/success_models_', combi, '.rds', sep = ''))
    return(list())
})

  
toc()


