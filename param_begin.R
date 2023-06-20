library(here)
library(fst)
library(ggplot2)
library(assertthat)
source(here("UtilityFunctions.R"))
library(tictoc)

library(tidyverse)
library(ggpubr)
library(readxl)
library(reshape2)
library(dbscan)
library(broom)
library(zoo)
library(knitr)
library(NMF,quietly = T)
library(MASS)
library(car)
library(doParallel)
library(foreach)
library(spatstat)
# library(raster)
library(gridExtra)
library(data.table)


# Original classification
# cells <- getCells()

# Alternative classification
# cells <- read_fst(here('DATA/SingleCells_altered.fst'))
# phenotype_list <- get_phenotypes()

# Use regional split

# structure_type <- 'Vascular stroma'
# structure_type <- 'Suppressed expansion'
# structure_type <- 'APC enriched'
# structure_type <- 'Granulocyte enriched'
structure_type <- "TLSlike"

print(structure_type)

cells <- readRDS(paste(here('DATA/splittedCells_'), gsub(' ','_', structure_type,fixed = T), '.rds',sep=''))


# Rename labels to fit methods of Alberto
cells <- cells %>% rename(tnumber = Split_ImageNumber)
cells <- cells %>% rename(phenotype = meta_description)
cells <- cells %>% rename(Xcenter = Location_Center_X)
cells <- cells %>% rename(Ycenter = Location_Center_Y)

# Remove unncessary attributes
cells <- cells[,c("tnumber", "ObjectNumber", "Xcenter", "Ycenter", "phenotype", "AreaShape_Area")]

tic('1-NN distances')
# call cores (parallelization)
# Comment this because of limitation on core usability in darwin
cores=detectCores()
cl <- makeCluster(floor(cores[1]/4))
registerDoParallel(cl)

# First NNdist all-against-all phenotypes
#' @param x Multiplex immunofluorescence coordinate dataframe (with Xcenter, Ycenter, phenotype, etc)
#' @param phenotype1 Phenotype From
#' @param phenotype2 Phenotype To
distances <- function(x, phenotype1,phenotype2) {
  # Create a point pattern
  spdat <- ppp(
    x = x$Xcenter,
    y = x$Ycenter,
    window = owin(
      xrange = c(min(x$Xcenter), max(x$Xcenter)),
      yrange = c(min(x$Ycenter), max(x$Ycenter))
    ),
    marks = x$phenotype
  )
  phenos <- unique(x$phenotype)
  
  #if phenotypes to compare are equal, then the 2nd nearest neighbour must be used, else the current point is compared to itself
  if(phenotype1 != phenotype2){
    res <- nncross(spdat[which(spdat$marks == phenotype1)], spdat[which(spdat$marks == phenotype2)], k = 1, what = 'dist')
  } else {
    res <- nncross(spdat[which(spdat$marks == phenotype1)], spdat[which(spdat$marks == phenotype2)], k = 2, what = 'dist')
    
  }
  return(res)
}

# Define which phenotypes will be studied (UNIQUE CELL TYPES IN YOUR DATASET)
phenotypes <- unique(cells$phenotype) #All cell types
# phenotypes <- unique(phenotype_list$stromal_immune$pg_id_map$meta_description)[1:2] #Take only 1 for debugging
print(paste0("The phenotypes we are considering: ", paste(phenotypes,collapse = ', ')))

# Compute 1-NN for all combinations of samples (n=24) x phenotype FROM (n=7) x phenotype TO (n=7)
spatresall <-
  foreach(i = unique(cells %>% pull(tnumber)),
          .packages = c('sp', 'spatstat'),
          .final = function(x) setNames(x, unique(cells %>% pull(tnumber)))) %dopar% {
            x <- cells[which(cells$tnumber == i),] 
            phenotypes <- setNames(phenotypes, phenotypes)
            spres <-
              lapply(phenotypes, function(p1) {
                lapply(phenotypes, function(p2) {
                  distances(x, p1, p2)
                })
              })
            return(spres)
          }

stopCluster(cl)

toc()


#Reshape hot mess of lists of lists of lists into long data frame.
dfl <- lapply(spatresall, function(tnr){
  lapply(tnr, function(p1){
    data.table::rbindlist(lapply(p1,function(Distance) {
      as.data.frame(Distance)
    }), idcol = "phenotype_to")
  })
})

dfl2 <- lapply(dfl, function(tnr){
  rbindlist(tnr, idcol = "phenotype_from")
})

dfl3 <- rbindlist(dfl2, idcol = "tnumber")

#add pseudocounts on distances for log transform
#1 pixel = 1 micron
dfl3$Distance <- (dfl3$Distance + 1)

NNdistances <- dfl3

saveRDS(NNdistances, file=here('scratch/NNdistances.rds'))

# free up some memory. these temps are no longer needed
rm(dfl)
rm(dfl2)
rm(dfl3)
gc() # do garbage collection to free up memory


options(dplyr.summarise.inform = FALSE)

tic('bin smoothing')
# Explore data from 0 to 300 microns
SlidingBins_300 <- rollapply(c(0:300),5,function(y) { return(c(min(y),max(y)))})

CountRows <- function(x) {
  cnt <- apply(SlidingBins_300,1, function(y) {nrow(x[which(x$Distance >= y[1] & x$Distance <= y[2])])})
  return(cnt)
}

BinCounts_300 <- apply(SlidingBins_300, 1, function(x) {
  CurrentWindow <- subset(NNdistances, Distance >= x[1] & Distance <= x[2])
  GroupCounts <- CurrentWindow %>% group_by(tnumber, phenotype_from, phenotype_to) %>%
    dplyr::summarise(N = dplyr::n())
  return(GroupCounts)
})

names(BinCounts_300) <- c(1:297)
dfslides_300 <- as.data.frame.matrix(SlidingBins_300)
dfslides_300$V3 <- apply(dfslides_300,1,mean)
dfslides_300$ID <- row.names(dfslides_300)

LongBinCounts_300 <- bind_rows(BinCounts_300, .id = "nth_window")
# WinMean is the X coordinate of your "calculated histogram"
LongBinCounts_300$WinMean <- dfslides_300$V3[match(LongBinCounts_300$nth_window, dfslides_300$ID)]

# # Area not available, so compute with max and min x and y coordinate
# # ALTERNATIVE METHOD? CONVEX HULL
# areas <- left_join(aggregate(Xcenter ~ tnumber, cells, function(x) min(x)),
#                    aggregate(Xcenter ~ tnumber, cells, function(x) max(x)), by='tnumber')
# areas <- left_join(areas,aggregate(Ycenter ~ tnumber, cells, function(x) min(x)), by='tnumber')
# areas <- left_join(areas,aggregate(Ycenter ~ tnumber, cells, function(x) max(x)), by='tnumber')
# areas$area <- (areas[,3]-areas[,2])*(areas[,5]-areas[,4])*0.001 #area in mm
# 
# # Match total areas of the tissue 
# LongBinCounts_300$tnumber <- as.integer(LongBinCounts_300$tnumber)
# LongBinCounts_300$Area <- LongBinCounts_300 %>%
#   left_join(data.frame(areas), by='tnumber') %>%
#   mutate(Area=area) %>% pull(Area)

# N.per.mm2 is the Y coordinate of your "calculated histogram"
# Normalize for area
# LongBinCounts_300$N.per.mm2 <- LongBinCounts_300$N / LongBinCounts_300$Area

# Might not be necessary due to AUC normalization [COMPARE]
LongBinCounts_300$N.per.mm2 <- LongBinCounts_300$N

# Function to compute Area under the curve 
funAUC <- function(x,y){
  id <- order(x)
  AUC <- sum(diff(x[id]) * rollmean(y[id],2))
  return(AUC)
}

AUCScaledSlides_300 <- LongBinCounts_300 %>%
  dplyr::group_by(phenotype_from, phenotype_to, tnumber) %>%
  dplyr::filter(n() >= 1) %>% 
  dplyr::mutate(N.per.mm2.scaled = N.per.mm2 / funAUC(WinMean, N.per.mm2))

print(paste('contains Inf values do to small number of cells:', (Inf %in% AUCScaledSlides_300$N.per.mm2.scaled)))

# Sometimes there are Infinite values here in case of small number of cells
# Remove these samples
AUCScaledSlides_300 <- AUCScaledSlides_300[!is.infinite(AUCScaledSlides_300$N.per.mm2.scaled),]

#double check if the AUC of AUCs is 1
AUCcheck <- AUCScaledSlides_300 %>%
  dplyr::group_by(phenotype_from, phenotype_to, tnumber) %>%
  dplyr::mutate(AUC = funAUC(WinMean, N.per.mm2.scaled))
stopifnot(length(unique(lapply(unique(AUCcheck$AUC), round)))==1)

#save this
write_delim(AUCScaledSlides_300 %>% 
              mutate(phenotype_combo=paste(phenotype_from,phenotype_to,sep='_to_')), paste(here('scratch/AUCScaledSlides_300_'),gsub(' ', '_', structure_type, fixed=T) ,'.tsv', sep=''),delim='\t')

toc()