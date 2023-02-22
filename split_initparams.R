library(here)

initial_parameters_firstrun_split <- split(initial_parameters_firstrun, f =  initial_parameters_firstrun$combo) 

for (c in names(initial_parameters_firstrun_split)){
  if (!file.exists(here(paste('scratch/init_parameters/init_params_', c, '.rds', sep='' )))){
    saveRDS(initial_parameters_firstrun_split[[c]], here(paste('scratch/init_parameters/init_params_', c, '.rds', sep='' )))
  }else{
    message(c, ' exists')
  }
}
