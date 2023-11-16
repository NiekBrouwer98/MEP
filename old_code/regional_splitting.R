library(tidyverse)
library(ggplot2)
library(here)
library(fst)
library(data.table)
source(here("UtilityFunctions.R"))
source(here('MEP_UtilityFunctions.R'))
library(ComplexHeatmap)
library(spatstat)
library(dbscan)

# structure_type <- 'Vascular stroma'
# structure_type <- 'Suppressed expansion'
# structure_type <- 'APC enriched'
# structure_type <- 'Granulocyte enriched'
# structure_type <- "TLSlike"
structure_type <- "FSP1+ enriched"
# structure_type <- "Active stroma"
# structure_type <- "PDPN+ active stroma"
# structure_type <- "CD8+ macrophages"
# structure_type <- "Active IR"

clinical_data <- getClinical()
cells <- getCellsAlternative()
structures <- getStructures()
structureLabels <- gsub("[^A-Za-z0-9+[:space:]]",'', structures$TME$labels)
structureLabels <- tibble(index = seq(1,10), label = structureLabels)
TMEStructures <- read_fst(here('scratch/TMEStructures.fst'), as.data.table=T)
allCommunities <- read_rds(here('scratch/allCommunities.rds'))

allCommunities <- merge(allCommunities, TMEStructures %>% dplyr::select(communityUID, TMEStructure), by='communityUID')

cellsWithoutExpression <- cbind(cells[,1:(ncol(cells)-43)],cells[,(ncol(cells)-3):ncol(cells)])
cellsWithStructure <- merge(cellsWithoutExpression, allCommunities %>% dplyr::select(-c('to', 'to_meta_id', 'nInteractions')), by.x = c('ImageNumber', 'ObjectNumber'), by.y=c('ImageNumber', 'from'), all=T)
cellsWithStructure <- cellsWithStructure %>% mutate(TMEStructure = as.numeric(TMEStructure))
cellsWithStructure <- merge(cellsWithStructure, structureLabels, by.x='TMEStructure', by.y='index',all.x=T) %>% distinct()

# Find image with a lot of vascular stroma
VS <- cellsWithStructure %>% filter(label == structure_type) %>% group_by(ImageNumber) %>% count(communityID) %>% count(ImageNumber)
VS_imageN <- VS %>% filter(n > 5) %>% pull(ImageNumber)
VS_images <- cellsWithStructure %>% filter(ImageNumber %in% VS_imageN)
VS_images <- VS_images %>% mutate(label =ifelse(label == structure_type, structure_type, 'other'))
VS_images$label[is.na(VS_images$label)] <- "other"

print(structure_type)

print(paste('number of images:', length(VS_imageN)))

# function to assign Tumor/Stroma using KDEtumor, KDEstrom, and input coordinates
structure_or_other <- function(kde_structure, kde_other, location){
  
  structure_value <- kde_structure %>% data.frame %>%
    mutate(dist=(x-location[1])^2+(y-location[2])^2) %>%
    arrange(dist) %>% slice(1) %>% pull(value)
  
  other_value <- kde_other %>% data.frame %>%
    mutate(dist=(x-location[1])^2+(y-location[2])^2) %>%
    arrange(dist) %>% slice(1) %>% pull(value)
  
  assigned_loc <- ifelse(structure_value > other_value, 'S', 'O')
  
  return(assigned_loc)
}

for(imageN in VS_imageN){
  # compute clusters
  dbclusters <- dbscan(as.matrix(VS_images %>% filter(ImageNumber == imageN) %>% dplyr::select(Location_Center_X,Location_Center_Y)), 
                       eps = 300, minPts = 50)
  # Assign clusters / foci's
  VS_images[VS_images$ImageNumber == imageN,'ClusterID'] <- dbclusters$cluster %>% as.character
}
VS_images <- VS_images %>% mutate(sample = paste(ImageNumber, ClusterID, sep='_'))

# DF to store optimal bandwiths:
sigmas_ppl <- data.frame(matrix(nrow=0,ncol=3))
colnames(sigmas_ppl) <- c('sample','sigma', 'structure_other')
sigmas_ppl$sample <- as.character(sigmas_ppl$sample)
sigmas_ppl$sigma <- as.numeric(sigmas_ppl$sigma)
sigmas_ppl$structure_other <- as.character(sigmas_ppl$structure_other)

samples <- unique(VS_images %>% pull(sample))
print(paste('Executing script for:',length(samples),'samples'))

for (i in samples){
  test <- VS_images %>% filter(sample == i) %>% dplyr::select(Location_Center_X, Location_Center_Y) %>% as.matrix

  # Compute point pattern for structure cells
  x <- VS_images %>% filter(sample == i) %>% filter(label == structure_type)
  spdat_structure <- ppp(
    x = x$Location_Center_X,
    y = x$Location_Center_Y,
    window = owin(
      xrange = c(min(x$Location_Center_X), max(x$Location_Center_X)),
      yrange = c(min(x$Location_Center_Y), max(x$Location_Center_Y))
    ),
    marks = x$label
  )
  # Window(spdat_structure) <- ripras(spdat_structure)

  # Compute point pattern for all other cells
  x <- VS_images %>% filter(sample == i) %>% filter(label != structure_type)
  spdat_other <- ppp(
    x = x$Location_Center_X,
    y = x$Location_Center_Y,
    window = owin(
      xrange = c(min(x$Location_Center_X), max(x$Location_Center_X)),
      yrange = c(min(x$Location_Center_Y), max(x$Location_Center_Y))
    ),
    marks = x$label
  )
  # Window(spdat_other) <- ripras(spdat_other)

  # Compute optimal sigma tumor
  optimal_sigma <- bw.ppl(spdat_structure)
  sigmas_ppl <- bind_rows(sigmas_ppl, 
                          data.frame(sample=i, sigma=(as.numeric(optimal_sigma)+(4*as.numeric(optimal_sigma))), structure_other='S'))
  print((as.numeric(optimal_sigma)+(4*as.numeric(optimal_sigma))))
  
  # Compute optimal sigma stroma
  optimal_sigma <- bw.ppl(spdat_other)
  sigmas_ppl <- bind_rows(sigmas_ppl, 
                          data.frame(sample=i, sigma=(as.numeric(optimal_sigma)), structure_other='O'))
  
  
  
}

colnames(sigmas_ppl) <- c('sample','sigma','structure_other')

overall_assigned_locs <- list()

for (i in samples){
  print(i)
  test <- VS_images %>% filter(sample == i) %>% dplyr::select(Location_Center_X, Location_Center_Y) %>% as.matrix
  ## Compute point pattern for tumor cells
  img <- VS_images %>% filter(sample == i) %>% filter(label == structure_type)
  spdat_structure <- ppp(
    x = img$Location_Center_X,
    y = img$Location_Center_Y,
    window = owin(
      xrange = c(min(test[,1]), max(test[,1])),
      yrange = c(min(test[,2]), max(test[,2]))
    ),
    marks = img$label
  )
  # Window(spdat_structure) <- ripras(spdat_structure)
  
  ## Compute point pattern for stromal cells
  img <- VS_images %>% filter(sample == i) %>% filter(label != structure_type)
  spdat_other <- ppp(
    x = img$Location_Center_X,
    y = img$Location_Center_Y,
    window = owin(
      xrange = c(min(test[,1]), max(test[,1])),
      yrange = c(min(test[,2]), max(test[,2]))
    ),
    marks = img$label
  )
  
  # Window(spdat_other) <- ripras(spdat_other)
  
  
  # Compute KDE structure
  sigma_structure <- sigmas_ppl %>% filter(sample == i) %>% 
    filter(get('structure_other') == 'S') %>% pull(sigma) %>% as.numeric
  k1_structure_raw <- density(spdat_structure, edge=TRUE, sigma=sigma_structure)
  # Compute KDE other
  sigma_other <- sigmas_ppl %>% filter(sample == i) %>% 
    filter(get('structure_other') == 'O') %>% pull(sigma) %>% as.numeric
  k1_other_raw <- density(spdat_other, edge=TRUE, sigma=sigma_other)
  
  # Classify Tumor / Stroma region for all cell positions from the sample
  ## KDE's are normalized by maximum value of KDE
  assigned_locs <- sapply(1:nrow(VS_images %>% filter(sample == i) ), 
                          function(x) structure_or_other(k1_structure_raw / max(k1_structure_raw), 
                                                         k1_other_raw / max(k1_other_raw), 
                                                         c(test[x,1], test[x,2])))
  overall_assigned_locs[[i]] <- assigned_locs
  
}


data_with_assigned_locations <- VS_images %>% slice(0) %>%
  mutate(ClusterID=as.character(ClusterID))

# Add assigned locations to each cell
for(i in samples){
  data_with_assigned_locations <- bind_rows(data_with_assigned_locations,
                                            VS_images %>% filter(sample == i) %>%
                                              mutate(assigned_loc = overall_assigned_locs[[i]])
  )
}

saveRDS(data_with_assigned_locations, paste(here('scratch/Regions_'),gsub(' ','_',gsub('+', 'plus',structure_type,fixed=T),fixed=T), '.rds', sep=''))

message('done!')
