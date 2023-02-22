#Generate a correlation matrix

library(ComplexHeatmap)

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

df_from <- all_parameters %>% select(c('phenotype_from', 'a', 'b'))
df_to <- all_parameters %>% select(c('phenotype_to', 'a', 'b'))

draw(Heatmap(df_from, name = "Pearson cor\n b coefficient", show_column_names=TRUE, 
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7, rot=45),
        cluster_rows  = TRUE, # turn off row clustering
        row_title_gp = gpar(font = c(1,1)),
        column_title="b parameter, CR",
        # top_annotation_height = unit(1, "mm"),
        show_row_dend = FALSE,show_column_dend =FALSE, column_names_side = "top",
        col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")) # scale color scale
)
)
