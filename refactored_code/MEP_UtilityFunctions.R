# This file contains all functions that are used over multiple notebooks
source(here("UtilityFunctions.R"))
library(tidyverse)

#'Print slide of given sample
#'
#'@return slide with highlighted cell type(s) and all cell types
show_slide <- function(sample_number, highlight_A){
  cells <- getCellsAlternative()
  cells <- cells %>% filter(meta_description %in% unique(cells$meta_description)[[1:10]])
  slide_1 <- ggplot(data=cells%>% filter(ImageNumber == sample_number) %>% mutate(type = ifelse(meta_description %in% highlight_A, 'highlighted', 'other'))) + geom_point(aes(x=Location_Center_X, y=Location_Center_Y, color=type), alpha =0.5) +
    theme_bw() + ggtitle(paste('imagenumber: ', sample_number))
  
  slide_2 <- ggplot(data=cells%>% filter(ImageNumber == sample_number)) + 
    geom_point(aes(x=Location_Center_X, y=Location_Center_Y, color=meta_description), alpha =0.5) + 
    theme_bw() + ggtitle(paste('imagenumber: ', sample_number)) + theme(legend.position = "none") + scale_color_manual(values = getDiscreteColors(10))
  
  return(plot_grid(slide_1, slide_2, n_col=2))
}

#'Print slide and distance distribution of given sample
#'
#'@return grid of two plots
show_distance_distribution <- function(cells, distance_data, sample, ylimit=0.2,color){
  sample_distance_data <- merge(x = distance_data %>% filter(tnumber %in% sample$tnumber) %>% filter(phenotype_combo %in%
                                                                                                            sample$phenotype_combo),
                                y = sample %>% dplyr::select(c(phenotype_combo, shape, scale)), all.X =T, by.x='phenotype_combo', by.y='phenotype_combo')
  
  dist <- ggplot(data= sample_distance_data %>% filter(tnumber == sample['tnumber'][[1]]) %>% filter(phenotype_combo == sample['phenotype_combo'][[1]])) + geom_bar(aes(x=distance_window, y=N.per.mm2.scaled), stat="identity",fill=color) +
    stat_function(fun = dweibull,colour= 'black', args = list(shape = sample['shape'][[1]], scale = sample['scale'][[1]]))+
    xlim(0,400) + ylim(0,ylimit) + theme_bw() + theme(legend.position = "none") + xlab('1-NN distance (micron)') + ylab('normalized counts')
  ''
  df <- cells%>% filter(ImageNumber == sample['tnumber'][[1]]) %>% 
    filter(meta_description == sample['phenotype_from'][[1]] | meta_description == sample['phenotype_to'][[1]]) %>% mutate(label = ifelse(meta_description == sample['phenotype_from'][[1]] , 'FROM', 'TO'))
  
  slide <- ggplot(data=df) +geom_point(aes(x=Location_Center_X, y=Location_Center_Y, color=label), alpha =1) +
    theme_bw() + DiscreteColors('label', fill=F, color=T)  + xlab('') + ylab('') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank())
  
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
  s <- readRDS(here('scratch/features/scale_parameters_run1.rds'))
  return(s)
}

getShapeFeatures <- function(){
  s <- readRDS(here('scratch/features/shape_parameters_run1.rds'))
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

getCellCounts <- function(cells){
  cell_counts <- cells %>% dplyr::select(c(ImageNumber, meta_description)) %>%
    dplyr::count(ImageNumber, meta_description)  %>%
    dplyr::rename(tnumber=ImageNumber)
  
  expand_cell_counts <- cell_counts %>% tidyr::expand(tnumber,meta_description)
  
  cell_counts <- merge(cell_counts, expand_cell_counts, by=c('tnumber', 'meta_description'), all.y=T) %>% 
    mutate(n = coalesce(n, 0))
  
  return(cell_counts)
}

getCombinationCounts <- function(cells){
  cell_counts_to <- getCellCounts(cells) %>% dplyr::rename(n_to = n)
  cell_counts_from <- getCellCounts(cells) %>% dplyr::rename(n_from = n)
  
  expand_cell_combinations <- cell_counts_to %>% tidyr::expand(tnumber, meta_description, meta_description) %>% 
    dplyr::rename(phenotype_from = meta_description...2) %>% 
    dplyr::rename(phenotype_to = meta_description...3)
  
  combination_counts <- merge(expand_cell_combinations, cell_counts_from, by.x=c('tnumber', 'phenotype_from'), by.y=c('tnumber', 'meta_description'))
  combination_counts <- merge(combination_counts, cell_counts_to, by.x= c('tnumber', 'phenotype_to'), by.y=c('tnumber', 'meta_description'))
  combination_counts <- combination_counts %>%
    unite('phenotype_combo', c(phenotype_from, phenotype_to), sep= '_to_', remove=F) %>%
    unite('unique_sample', c(tnumber, phenotype_combo), remove=F)
  
  return(combination_counts)
}

getTumorAndTMETypes <- function(){
  cells <- getCells()
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

getNonTumourImages <- function(){
  cells <- getCells()
  return(unique(cells %>% filter(isTumour == F) %>% pull(ImageNumber)))
}

getCellsAlternative <- function(){
  cells <- read_fst(here('DATA/SingleCells_altered.fst'), as.data.table = T)
  return(cells)
}

getCombinationsAlternative <- function(){
  combinations <- readRDS(here('scratch/combinations_selection_run2.rds'))
  return(combinations)
}

getALLParametersAlternative <- function(){
  all_parameters <- read_rds(here('scratch/all_parameters_run2.rds'))
  return(all_parameters)
  
}

getScaleFeaturesAlternative <- function(){
  s <- readRDS(here('scratch/features/scale_parameters_run3.rds'))
  return(s)
}

getShapeFeaturesAlternative <- function(){
  s <- readRDS(here('scratch/features/shape_parameters_run3.rds'))
  return(s)
}

getDensityFeaturesAlternative <- function(){
  s <- readRDS(here('scratch/features/cell_counts_per_imageAlternative.rds'))
  return(s)
}

getCellProportionsPerImageAlternative <- function(){
  s <- readRDS(here('scratch/features/cell_proportions_per_imageAlternative.rds'))
  return(s)
}
