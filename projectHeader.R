# Header file

# libraries
loadLibrary <- function(x) { 
  suppressMessages(library(x, character.only = T, quietly = T))
}
libraries <- c("data.table", "ggplot2", "parallel", 
  "assertthat", "RColorBrewer", "fst", 'latex2exp', 'here')
invisible(lapply(libraries, loadLibrary))

# project specific functions are housed in -./code/projectPkg- library;
# build and load using roxygen2
loadprojectPkg <- function(){
  if('projectPkg' %in% (.packages()) == FALSE) {
    roxygen2::roxygenise(package.dir = here('Code/projectPkg'))
  } else{
    cat('projectPkg loaded\n')
  }
}
suppressWarnings(loadprojectPkg())
rm(libraries, loadprojectPkg, loadLibrary)

# Load panel
projectPanel <- getPanel()

# set seed
projectSeed <- 89230689
set.seed(projectSeed)

# set working directory as outdir
if(!dir.exists(here('scratch'))) dir.create(here('scratch'))
setwd(here('scratch'))

# nThreads
n_threads <- parallel::detectCores()
setDTthreads(n_threads)

# ids vars
idvars <- c('ImageNumber', 'ObjectNumber')

# display width
options(width = 180)