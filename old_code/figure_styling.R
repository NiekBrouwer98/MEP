# This file contains all the functions to style figures for reporting
library(rcartocolor)
library(palettes)
library(circlize)
source(here("UtilityFunctions.R"))
source(here('MEP_UtilityFunctions.R'))

# "#88CCEE" "#CC6677" "#DDCC77" "#888888"

alternativeCellTypeLabels <- function(){
  cellType_labels <- tibble(name = unique(getCellsAlternative() %>% pull(meta_description))) %>% arrange(name)
  cellType_labels <- cellType_labels %>% mutate(conversion = paste('$',name, '$',sep=''))
  # cellType_labels <- cellType_labels %>% mutate(conversion = gsub("&","$&$", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("ERorHER2-","", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("Basal","Basal_cells", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("Endothelial"," Endothelial", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("Fibroblasts"," Fibroblasts", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("Myofibroblasts"," Myofibroblasts ", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("B cells"," B cells", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("Ki67+_cells_epithelial","epithelial_Ki67+_cells", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("Ki67+_cells_TME","TME_Ki67+_cells", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub(" ","$ $", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("+","^{+}", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("-","^{-}", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("-","^{-}", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("med","^{med}", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("_","$ $", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("Granulocytes$ $Macrophages","Granulocytes$ $$&$$ $Macrophages", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("_A","$&$ A", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("_","$ $", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("or","", conversion, fixed=T))

  
  return(cellType_labels)
}


originalCellTypeLabels <- function(){
  cellType_labels <- tibble(name = unique(getCells() %>% pull(meta_description))) %>% arrange(name)
  cellType_labels <- cellType_labels %>% mutate(conversion = paste('$',name, '$',sep=''))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("&","$&$", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub(" ","$ $", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("-","^{-}", conversion, fixed=T))
  cellType_labels <- cellType_labels %>% mutate(conversion = gsub("CK8^{-}18","CK8-18", conversion, fixed=T))
  
  
  return(cellType_labels)
  
}


heatmapColorScale <- function(min, max){
  # col_fun = colorRamp2(c(min, 0, max), c( viridis::viridis(3)[[1]],viridis::viridis(3)[[2]],viridis::viridis(3)[[3]]))
  #  carto_pal(3, "Tropic")
  col_fun = colorRamp2(c(min, 0, max), c(carto_pal(4, "Safe")[[1]], 'white', carto_pal(4, "Safe")[[2]]))
  
  return(col_fun)
}

heatmapfraction <- function(min,max){
  col_fun = colorRamp2(c(min, max), c('white',carto_pal(4, "Safe")[[2]]))
  return(col_fun)

}

DiscreteColors <- function(legend_name,fill=T, color=F){
  if(fill){
    return(scale_fill_carto_d(name = legend_name, palette = "Safe"))
  }
  
  if(color){
    return(scale_color_carto_d(name = legend_name, palette = "Safe"))
  }
  
  else{
    return(NULL)
  }

}

getDiscreteColors <- function(n){
  return(carto_pal(n, "Safe"))
}

getDivergingColors <- function(n){
  return(carto_pal(n,'Geyser'))
}

getSequentialColors <- function(n){
  print('pick odd number')
  return(carto_pal(n,'Geyser')[ceiling(n/2):(n)])
}
