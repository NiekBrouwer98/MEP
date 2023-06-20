# This file contains all the functions to style figures for reporting
library(rcartocolor)
library(palettes)
library(circlize)
source(here("UtilityFunctions.R"))
source(here('MEP_UtilityFunctions.R'))

alternativeCellTypeLabels <- function(){
  cellType_labels <- tibble(name = unique(getCellsAlternative() %>% pull(meta_description))) %>% arrange(name)
  cellType_labels <- cellType_labels %>% mutate(conversion = paste('$',name, '$',sep=''))
  # cellType_labels <- cellType_labels %>% mutate(conversion = gsub("&","$&$", conversion, fixed=T))
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
  col_fun = colorRamp2(c(min, 0, max), c( 'blue','white','red'))
  
  return(col_fun)
}

heatmapfraction <- function(min,max){
  col_fun = colorRamp2(c(min, max), c(viridis::viridis(2)[[1]], viridis::viridis(2)[[2]]))
  return(col_fun)

}

categoricalColorScale <- function(legend_name,fill=T, color=F,color_type){
  if(fill){
    return(scale_fill_carto_c(name = legend_name, palette = "ag_Sunset"), type=color_type)
  }
  
  if(color){
    return(scale_color_carto_c(name = legend_name, palette = "ag_Sunset"), type=color_type)
  }
  
  else{
    return(NULL)
  }
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