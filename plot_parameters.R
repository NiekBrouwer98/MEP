# Make plots of the parameters

library(here)
library(ggplot2)
library(tidyverse)

combinations <- readRDS(here('scratch/combinations_selection.rds'))
toCD8_combinations <-  grep("to_CD8^{+} T cells", combinations, value = TRUE, fixed= TRUE)
toCD4APC_combinations <- grep("to_CD4^{+} T cells & ", combinations, value = TRUE, fixed=TRUE)
toCD4_combinations <- setdiff(grep("to_CD4^{+} T cells", combinations, value = TRUE, fixed=TRUE), toCD4APC_combinations)
toFoxP3_combinations <-  grep("to_T_{Reg}", combinations, value = TRUE, fixed= TRUE)
toMacrophage_combinations <-  grep("to_Macrophages$", combinations, value = TRUE)
toBcell_combinations <-  grep("to_B cells", combinations, value = TRUE, fixed= TRUE)
# toHER2_combinations <- grep("to_HER2^{+}$", combinations, value = TRUE, fixed= TRUE)
# toCK818ER_combinations <- grep("to_CK8-18^{+} ER^{hi}$", combinations, value = TRUE, fixed= TRUE)
toCK818CXCL12_combinations <- grep("to_CK8-18^{hi}CXCL12^{hi}", combinations, value = TRUE, fixed =TRUE) #Cancer

get_parameters <- function(combinations){
  params <- c()
  for (c in combinations){
    tryCatch({
    success_model <- readRDS(here(paste('scratch/success_models/success_models_', c,'.rds', sep='')))
    params <- append(params, list(success_model[[3]]))
    names(params)[length(params)] <- c
    },error=function(e){print(c)})
    
  }
  
  params_result <- bind_rows(params, .id='phenotype_combo')
  return(params_result)
}

all_parameters <- get_parameters(combinations) %>%
  separate(phenotype_combo, into=c('phenotype_from','phenotype_to'),sep='_to_') %>%
  filter(a > 0.5) %>%
  filter(b > 8)

means <- all_parameters %>% 
  group_by(phenotype_from, phenotype_to ) %>% 
  summarise(mean_a = mean(a),
            mean_b  = mean(b))
  

library(ggrepel)

ggplot(data = means, aes(x=mean_a, y=mean_b, label=paste(phenotype_from, '_to_', phenotype_to, sep=''))) +
  geom_point(data = all_parameters, aes(x=a, y=b), alpha=0.5, size=0.5) +
  geom_point(data = means, size = 2, aes(x=mean_a, y=mean_b, label=paste(phenotype_from, '_to_', phenotype_to, sep='')), color='red') +
  # scale_y_log10() +
  # geom_text_repel(data          = subset(means,mean_a > quantile(mean_a, 0.99)),
  #                 nudge_y       = 600,
  #                 size          = 5,
  #                 box.padding   = 1.5,
  #                 point.padding = 0.5,
  #                 force         = 200,
  #                 segment.size  = 1,
  #                 segment.color = "grey50",
  #                 direction     = "both") +
  geom_text_repel(data          = subset(means,mean_b > quantile(mean_b, 0.99)),
                  nudge_y       = 600,
                  # nudge_x       = 7.5,
                  size          = 5,
                  box.padding   = 1.5,
                  point.padding = 0.5,
                  force         = 200,
                  segment.size  = 1,
                  segment.color = "grey50",
                  direction     = "both") +
  xlab('Shape') + ylab('Scale')


ggsave(paste(here('output/'), 'allparameters_withoutliersb.png', sep=''))

plot_parameters <- function(combi){
  point <- all_parameters %>%
    filter(phenotype_to == combi) %>%
    ggplot(aes(x=a, y=b, color=phenotype_from)) +
    geom_point(alpha=0.7, size=3) +
    ylim(10,300) + xlim(0,6) +
    theme_bw() + scale_y_log10() +
    xlab('Shape') + ylab('Scale') +
    theme(legend.position='right') + ggtitle(paste('From X to ', combi)) + guides(colour = guide_legend(override.aes = list(size=10)))
  
  print(point)
  
  ggsave(paste(here('output/'), combi, '.png', sep=''))
}

six_combinations <- c('CK8-18^{hi}CXCL12^{hi}', 'B cells','CD4^{+} T cells',
                      'Macrophages','CD8^{+} T cells','T_{Reg} & T_{Ex}')

for (c in six_combinations){
  plot_parameters(c)
}


plot_sixcelltypes <- function(combi){
  point <- all_parameters %>%
  filter(phenotype_to == combi) %>%
  filter(phenotype_from == 'CD8^{+} T cells' | phenotype_from == 'CD4^{+} T cells' |
         phenotype_from == 'Macrophages' | phenotype_from == 'B cells' |
           phenotype_from == 'CK8-18^{hi}CXCL12^{hi}' | phenotype_from == 'T_{Reg} & T_{Ex}' ) %>%
  ggplot(aes(x=a, y=b, color=phenotype_from)) +
  geom_point(alpha=0.7, size=2) +
  ylim(10,300) + xlim(0,6) +
  theme_bw() + scale_y_log10() +
  xlab('Shape') + ylab('Scale') + scale_color_brewer(palette='Set2') +
  theme(legend.position='right') + ggtitle(paste('From X to ', combi))
  
  print(point)

  
  ggsave(paste(here('output/'), combi, '_sixcelltypes.png', sep=''))
}

for (c in six_combinations){
  plot_sixcelltypes(c)
}

