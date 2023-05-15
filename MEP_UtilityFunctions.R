# This file contains all functions that are used over multiple notebooks
source(here("UtilityFunctions.R"))
library(tidyverse)
cells <- getCells()

#'Print slide of given sample
#'
#'@return slide with highlighted cell type(s) and all cell types
show_slide <- function(sample_number, highlight_A){
  slide_1 <- ggplot(data=cells%>% filter(ImageNumber == sample_number) %>% mutate(type = ifelse(meta_description %in% highlight_A, 'highlighted', 'other'))) + geom_point(aes(x=Location_Center_X, y=Location_Center_Y, color=type), alpha =0.5) +
    theme_bw() + ggtitle(paste('imagenumber: ', sample_number))
  
  slide_2 <- ggplot(data=cells%>% filter(ImageNumber == sample_number)) + 
    geom_point(aes(x=Location_Center_X, y=Location_Center_Y, color=meta_description), alpha =0.5) + 
    theme_bw() + ggtitle(paste('imagenumber: ', sample_number))
  
  return(plot_grid(slide_1, slide_2, n_col=2))
}

#'Print slide and distance distribution of given sample
#'
#'@return grid of two plots
show_distance_distribution <- function(sample, ylimit=0.02){
  sample_distance_data <- merge(x = all_distances_data %>% filter(tnumber %in% sample$tnumber) %>% filter(phenotype_combo %in%
                                                                                                            sample$phenotype_combo),
                                y = sample %>% select(c(phenotype_combo, shape, scale)), all.X =T, by.x='phenotype_combo', by.y='phenotype_combo')
  
  dist <- ggplot(data= sample_distance_data %>% filter(tnumber == sample['tnumber'][[1]]) %>% filter(phenotype_combo == sample['phenotype_combo'][[1]])) + geom_bar(aes(x=distance_window, y=N.per.mm2.scaled), stat="identity",alpha=0.5) +
    stat_function(fun = dweibull, args = list(shape = sample['shape'][[1]], scale = sample['scale'][[1]]))+
    xlim(0,300) + ylim(0,ylimit) + theme_bw() + theme(legend.position = "none") + ggtitle(label = sample['phenotype_combo'][[1]])
  
  slide <- ggplot(data=cells%>% filter(ImageNumber == sample['tnumber'][[1]]) %>% 
                    filter(meta_description == sample['phenotype_from'][[1]] | meta_description == sample['phenotype_to'][[1]])) + 
    geom_point(aes(x=Location_Center_X, y=Location_Center_Y, color=meta_description), alpha =0.5) +
    theme_bw() + ggtitle(paste('imagenumber: ', sample['tnumber'][[1]]))
  
  return(plot_grid(dist, slide,ncol=2))
  
}

#'Generate matrix from dataframe with filtering options
#'
#'@return matrix as dataframe
generate_matrix <- function(df, col_rownames ,subselection=NULL, NA_percentage=0, scale_rows=FALSE){
  m <- as.data.frame(df)
  rownames(m) <- m[,col_rownames]
  m <- subset(m,select = -c(get(col_rownames)))
  
  # Filtering
  m <- m[, which(colMeans(!is.na(m)) > NA_percentage)]
  
  if (!is.null(subselection)){
    m <- m %>% select(all_of(subselection))
  }
  
  # Scaling
  m <- scale(m)
  
  if (scale_rows == TRUE){
    m <- t(scale(t(m)))
  }
  m <- replace(m,is.na(m),0)
  
  m <- as.data.frame(m)
  
  return(m)
}

#' The first parameter estimation was done with thresholds [0,20,50,70,100]
#' The second parameter estimation was done with thresholds [0,2,5,10,20]
#' Merge the two estimations
#' @return Save parameter file in scratch dir
CompileALLParameters <- function(){
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
  
  get_all_parameters <- function(save){
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
  
  all_parameters <- get_all_parameters(save=TRUE)
  
}



transform_parameters_to_matrix <- function(outliers = TRUE){
  if (outliers == TRUE){
    all_parameters <- getALLParameters()
    
  }else{
    all_parameters <- readRDS(here('scratch/all_parameters_withoutOutliers.rds'))
  }
  shape_parameters <- tibble(tnumber = unique(all_parameters$tnumber))
  scale_parameters <- tibble(tnumber = unique(all_parameters$tnumber))
  
  for (c in unique(all_parameters$phenotype_combo)){
    shape_parameters <- left_join(x = shape_parameters, y= all_parameters %>% 
                                    select(c(tnumber, phenotype_combo, shape)) %>% 
                                    filter(phenotype_combo == c),
                                  by='tnumber') %>%
      select(-c(phenotype_combo))
    names(shape_parameters)[names(shape_parameters) == 'shape'] = c
    
    scale_parameters <- left_join(x = scale_parameters, y= all_parameters %>% 
                                    select(c(tnumber, phenotype_combo, scale)) %>% 
                                    filter(phenotype_combo == c),
                                  by='tnumber') %>%
      select(-c(phenotype_combo))
    names(scale_parameters)[names(scale_parameters) == 'scale'] = c
  }
  
  saveRDS(shape_parameters, here('scratch/features/shape_parameters.rds'))
  saveRDS(scale_parameters, here('scratch/features/scale_parameters.rds'))
}


addLabelToGrid <- function(plotGrid, title){
  title <- ggdraw() + 
    draw_label(title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  plotWithTitle <- plot_grid(
    title, plotGrid,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  return(plotWithTitle)
}

### GET functions

getStructureFeatures <- function(){
  s <- readRDS(here('scratch/features/structure_proportions_per_image.rds'))
  return(s)
}

getNetworkFeatures <- function(){
  s <- readRDS(here('scratch/features/network_properties_per_patient.rds'))
  return(s)
}

getScaleFeatures <- function(){
  s <- readRDS(here('scratch/features/scale_parameters.rds'))
  return(s)
}

getShapeFeatures <- function(){
  s <- readRDS(here('scratch/features/shape_parameters.rds'))
  return(s)
}

getDensityFeatures <- function(){
  s <- readRDS(here('scratch/features/cell_counts_per_image.rds'))
  return(s)
}

getCellProportionsPerPatient <- function(){
  s <- readRDS(here('scratch/features/cell_proportions_per_patient.rds'))
  return(s)
}

getCellProportionsPerImage <- function(){
  s <- readRDS(here('scratch/features/cell_proportions_per_image.rds'))
  return(s)
}

getCellCounts <- function(){
  cell_counts <- cells %>% select(c(ImageNumber, meta_description)) %>%
    dplyr::count(ImageNumber, meta_description)  %>%
    rename(tnumber=ImageNumber)
  
  expand_cell_counts <- cell_counts %>% tidyr::expand(tnumber,meta_description)
  
  cell_counts <- merge(cell_counts, expand_cell_counts, by=c('tnumber', 'meta_description'), all.y=T) %>% 
    mutate(n = coalesce(n, 0))
  
  return(cell_counts)
}

getCombinationCounts <- function(){
  cell_counts_to <- getCellCounts() %>% rename(n_to = n)
  cell_counts_from <- getCellCounts() %>% rename(n_from = n)
  
  expand_cell_combinations <- cell_counts_to %>% tidyr::expand(tnumber, meta_description, meta_description) %>% 
    rename(phenotype_from = meta_description...2) %>% 
    rename(phenotype_to = meta_description...3)
  
  combination_counts <- merge(expand_cell_combinations, cell_counts_from, by.x=c('tnumber', 'phenotype_from'), by.y=c('tnumber', 'meta_description'))
  combination_counts <- merge(combination_counts, cell_counts_to, by.x= c('tnumber', 'phenotype_to'), by.y=c('tnumber', 'meta_description'))
  combination_counts <- combination_counts %>%
    unite('phenotype_combo', c(phenotype_from, phenotype_to), sep= '_to_', remove=F) %>%
    unite('unique_sample', c(tnumber, phenotype_combo), remove=F)
  
  return(combination_counts)
}

getTumorAndTMETypes <- function(){
  TME =  c("B cells","CD38^{+} lymphocytes","CD4^{+} T cells","CD4^{+} T cells & APCs","CD57^{+}","CD8^{+} T cells","Endothelial","Fibroblasts","Fibroblasts FSP1^{+}","Granulocytes", "Ki67^{+}","Macrophages","Macrophages & granulocytes", "Myofibroblasts","Myofibroblasts PDPN^{+}", "T_{Reg} & T_{Ex}")
  tumor = setdiff(unique(cells %>% pull(meta_description)), TME)
  
  return(list(TME_types=TME, tumor_types=tumor))
}

getCombinations <- function(){
  combinations <- readRDS(here('scratch/combinations_selection.rds'))
  return(combinations)
}

getALLParameters <- function(){
  all_parameters <- read_rds(here('scratch/all_parameters.rds'))
  return(all_parameters)
}

getOutliers <- function(){
  outliers <- read_rds(here('scratch/outlier_parameters.rds'))
  reutrn(outliers)
}
