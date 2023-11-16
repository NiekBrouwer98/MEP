# Compute median distances for every cell type-cell type combination
# Input: 'scratch/NNdistances.rds' containing the distance distributions computed in alberto_implementation.Rmd

library(here)
library(foreach)

all_distances_data <- read_rds(here('scratch/AUCScaledSlides_300_ALL.rds'))
combinations <- unique(all_distances_data %>% pull(phenotype_combo))


median_distances <- foreach(c= 1:length(combinations), .combine = 'rbind') %dopar% {
  samples <- unique(all_distances_data %>% filter(phenotype_combo == combinations[c]) %>% pull(tnumber))
  for (sample in samples){
    
  }
  distances <- all_distances_data %>% filter(tnumber %in% sample$tnumber) %>% filter(phenotype_combo %in%
                                                                          sample$phenotype_combo)
}

