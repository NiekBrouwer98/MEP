library(ggplot2)
library(tidyverse)
library(gridExtra)

all_distances_data <- readRDS(file  = here('scratch/AUCScaledSlides_300_ALL.rds'))
initial_params_result <- readRDS(file = here('scratch/initial_parameters.rds'))

all_distances_data <- all_distances_data %>% as_tibble %>%
  mutate(distance_window = WinMean) %>%
  mutate(phenotype_combo = paste(phenotype_from, phenotype_to, sep='_to_')) %>%
  dplyr::select(tnumber, phenotype_combo, `N.per.mm2.scaled`, distance_window) %>%
  mutate(new = `N.per.mm2.scaled` * 1000) %>%
  mutate(new =round(new))

# combinations <- unique(all_distances_data %>% pull(phenotype_combo))
combinations_selection <- readRDS(file = here('scratch/combinations_selection.rds'))




plot_histogram = function(t){
  combinations_t <- unique(filter(initial_params_result, tnumber==t)['combo'][[1]])
  combinations_sample <- sample(1:length(combinations_t), 1, replace=FALSE)
  c = combinations_t[[combinations_sample]]
  return(ggplot(all_distances_data[all_distances_data$tnumber == t & all_distances_data$phenotype_combo == c ,]) +
    geom_bar(aes(x=distance_window, y=N.per.mm2.scaled), stat="identity") +
    stat_function(fun = dweibull, args = list(shape = filter(initial_params_result, tnumber==t & combo==c & term=='shape')['estimate'][[1]],
                                              scale = filter(initial_params_result, tnumber==t & combo==c & term=='scale')['estimate'][[1]])) + ggtitle(paste("tnumber: ", t, "combi: ", c)) + xlim(0,300))
}

plts <- lapply(sample(1:length(unique(all_distances_data$tnumber)), 3, replace=FALSE), plot_histogram)
do.call(grid.arrange, plts)


plot_notexistingparams = function(t){
  combinations_t <- setdiff(combinations_selection, unique(filter(initial_params_result, tnumber==t)['combo'][[1]]))
  combinations_sample <- sample(1:length(combinations_t), 1, replace=FALSE)
  c = combinations_t[[combinations_sample]]
  return(ggplot(all_distances_data[all_distances_data$tnumber == t & all_distances_data$phenotype_combo == c ,]) +
           geom_bar(aes(x=distance_window, y=N.per.mm2.scaled), stat="identity") +
           stat_function(fun = dweibull, args = list(shape = filter(initial_params_result, tnumber==t & combo==c & term=='shape')['estimate'][[1]],
                                                     scale = filter(initial_params_result, tnumber==t & combo==c & term=='scale')['estimate'][[1]])) + ggtitle(paste("tnumber: ", t, "combi: ", c)) + xlim(0,300))
}

plts_notexistingparams <- lapply(sample(1:length(unique(all_distances_data$tnumber)), 3, replace=FALSE), plot_notexistingparams)
do.call(grid.arrange, plts_notexistingparams)

plot_specific = function(t){
  c = 'CK8-18^{+} ER^{hi}_to_CK8-18^{+} ER^{hi}'
  return(ggplot(all_distances_data[all_distances_data$tnumber == t & all_distances_data$phenotype_combo == c ,]) +
           geom_bar(aes(x=distance_window, y=N.per.mm2.scaled), stat="identity") +
           stat_function(fun = dweibull, args = list(shape = filter(initial_params_result, tnumber==t & combo==c & term=='shape')['estimate'][[1]],
                                                     scale = filter(initial_params_result, tnumber==t & combo==c & term=='scale')['estimate'][[1]])) + ggtitle(paste("tnumber: ", t, "combi: ", c)) + xlim(0,300))
}

plts_notexistingparams <- lapply(sample(1:length(unique(all_distances_data$tnumber)), 3, replace=FALSE), plot_specific)
do.call(grid.arrange, plts_notexistingparams)
