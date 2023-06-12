library(fst)

#' Load antibody panel details
#' 
#' Loads and formats antibody panel
#'@import data.table 
#'@return List of panel table and character vectors of channel
#'  names (epithelial, stromal_immune and all_channels) 
getPanel <- function(){
  
  # Utility function for renaming
  nicelabs <- function(t, name, label, clone = NULL) {
    if (!is.null(clone)) {
      s_row <-
        intersect(grep(t, panel[, target]), grep(clone, panel[, antibodyclone]))
    } else{
      s_row <- grep(t, panel[, target])
    }
    panel[s_row, var_name := name]
    return(panel[s_row, var_label := label])
  }
  
  # Get data and generate vars
  panel <- here('AbPanel.csv')
  panel <- fread(panel)
  panel[, var_name := tolower(target)]  
  nicelabs("Histone H3", "hh3", "Histone H3")  
  nicelabs("Cytokeratin 5", "ck5", "CK5")
  nicelabs("Cytokeratin 8/18", "ck8_18", "CK8-18")
  nicelabs("CD278","icos","ICOS")
  nicelabs("CD134","ox40","OX40")
  nicelabs("CD279","pd1", "PD-1")
  nicelabs("GITR","gitr","GITR")
  nicelabs("Rabbit","er","ER")
  nicelabs("CD140b","pdgfrb","PDGFRB")
  nicelabs("CD31_vWF","cd31_vwf","CD31-vWF")
  nicelabs("CXCL12","cxcl12","CXCL12")
  nicelabs("Beta-2","b2m","B2M")
  nicelabs("pan","pan_ck","panCK")
  nicelabs("Cleaved","c_c_3","c-Caspase3")
  nicelabs("CD8a", "cd8", "CD8")

  nicelabs("c-erbB-2", "her2_3b5", "HER2 (3B5)", "3B5")  
  nicelabs("c-erbB-2", "her2_d8f12", "HER2 (D8F12)", "D8F12")  
  
  panel[metaltatC == "Ir191", var_name := "DNA1"]  
  panel[metaltag == "Ir193", var_name := "DNA2"]
  panel[is.na(var_label), var_label := target]
  panel[var_label=="", var_label := var_name]
  
  panel[target == "CD45", var_name := "cd45ro"]
  panel[target == "CD45", var_label := "CD45RO"]
  
  panel[, var_name := gsub("-", "_", var_name)]
  epithelial <- c("h3", "her2_3b5", "her2_d8f12", "b2m", "cd57", 
                  "ki_67", "ck5", "hla_abc", "pan_ck", "c_c_3", 
                  "hla_dr", "ck8_18", "cd15", "er", "cxcl12", 'sma', 
                  'podoplanin') # sma for myoepithelial
  stromal_immune <- c("hh3", "icos", "ox40", "cd68", "cd3", "podoplanin", "cd11c", "pd1", 
    "gitr", "cd16", "sma", "cd45ra", "b2m", "cd45ro", "foxp3", "cd20", "cd8", "cd57", "ki_67", "caveolin_1", 
    "cd4", "cd31_vwf", "hla_abc", "c_c_3", "cd38", "hla_dr", "cd15", "fsp1", "cd163", "pdgfrb") 
  
  panel[, epithelial := grepl(paste0(epithelial, collapse = "|"), var_name) ]
  panel[, stromal_immune := grepl(paste0(stromal_immune, collapse = "|"), var_name) ]
  functionOrder <- panel[order(function_order), function_order]
  names(functionOrder) <- panel[order(function_order), var_label]
  functionOrder <- functionOrder[!is.na(functionOrder)]

  # Output
  setkey(panel, to_sort)
  stromal_immune <- panel[stromal_immune == T, var_label]
  epithelial <- panel[epithelial == T, var_label]
  projectPanel <- list(
    panel_dat = panel[, .(metaltag, var_name, var_label, epithelial, stromal_immune)],
    all_channels = panel[, var_label],
    stromal_immune = stromal_immune,
    epithelial = epithelial,
    function_order = functionOrder
  )

  return(projectPanel)

}  


#'Load formatted clinical data
#'
#'@return data.table of formatted clincial data
getClinical <- function() {
  clinical <- read_csv(here('DATA/ClinicalMerged.csv'))
  cells <- getCells()
  imagenumber_metaid <- cells %>% distinct(ImageNumber, metabric_id)
  
  clinical_data <- merge(clinical, imagenumber_metaid, by='metabric_id', all.y=T)
  
  return(clinical_data)
}


#'Load formatted single cell data
#'
#'@return data.table of formatted single cell data
getCells <- function() {
  cells <- read_fst(here('DATA/SingleCells.fst'), as.data.table = T)
  # cells <- cells %>% mutate(meta_description = gsub("[^A-Za-z0-9+\\-]", "", meta_description))
  return(cells)
}


#' Get tumour subtype colours
#' 
#' @return named list (PAM50 and IntClust) of colours 
getSubtypeCols <- function(){
  pam50 <-
    c(
    'Luminal A' = "#1F78B4",  
    'Luminal B' = "#A6CEE3",    
    'HER2' = "#FB9A99",  
    'Basal' = "#E41A1C",  
    'Normal-like' = "#66A61E"
    )  

  intclust <- 
    c(
    'IntClust 1' = '#FF5500',
    'IntClust 2' = '#00EE76',
    'IntClust 3' = '#CD3278',
    'IntClust 4+' = '#00C5CD',
    'IntClust 4-' = '#B5D0D2',
    'IntClust 5+' = '#8B0000',
    'IntClust 5-' = '#D08282', 
    'IntClust 6' = '#FFFF40',
    'IntClust 7' = '#0000CD',
    'IntClust 8' = '#FFAA00',
    'IntClust 9' = '#EE82EE',
    'IntClust 10' = '#7D26CD'
    )
  
  molecular <- 
    c(
      'ER+HER2+' = '#FF5500',
      'ER-HER2+' = '#00EE76',
      'ER+HER2-' = '#CD3278',
      'ER-HER2-' = '#00C5CD'

    )

  return(list(
    PAM50 = pam50,
    IntClust = intclust,
    molecularSubtypes = molecular
  ))
}


#' Project ggplot theme
#' Thanks to https://github.com/taylor-lab/composite-mutations/blob/master/r/prerequisites.R
#' @param base_size defaults to 11; cascades to line and rect unless otherwise set
#' @param base_line_size defaults to base_size/22
#' @param base_rect_size defaults to base_size/22
theme_prj <- function(base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
          line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"),
          text = element_text(family = 'ArialMT', face = "plain",
                              colour = "black", size = base_size, lineheight = 0.9,
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F),
          axis.text = element_text(colour = "black", family='ArialMT', size=rel(0.8)),
          axis.ticks = element_line(colour = "black", size=rel(1)),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = rel(1)),
          legend.key = element_blank(),
          strip.background = element_blank()
    )
}


#' Convert ggplot into two objects - plot and legend
#' Thanks to https://github.com/taylor-lab/composite-mutations/blob/master/r/prerequisites.R
#' @param p a ggplot2 object
#' @return named list of 2 - plot and legend
#' @import cowplot ggplot2 grid gtable
extract_gglegend <- function(p){

    ## extract the legend from a ggplot object
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if(length(leg) > 0) leg <- tmp$grobs[[leg]]
    else leg <- NULL
    leg

    ## return the legend as a ggplot object
    legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
    plot <- p + theme(legend.position='none')
    list(plot=plot,legend=legend)
}


#' Utility to parse latex expressions as axis labels
#' To use e.g. \code{scale_y_discrete(labels = label_parse)}
#' @import ggplot2 latex2exp
label_parse <- function(breaks) {parse(text = latex2exp::TeX(breaks))}


#' Scale (z-score) and clip
#' 
#' Z-scores a vector and clips min and max at \code{limit} 
#' @param x numeric vector
#' @param limit scalar value to clip at; default is 2
#' @return numeric vector
scale_clip <- function(x, limit = 2){
  x <- scale(x)
  x[x > limit] <- limit
  x[x < -limit] <- -limit
  return(x)
}


#' Cell phenotype characteristics
#'
#' Utility function to return meta phenotypes, meta ids and colours
#' @param TeX_to_expression boolean, defaults to false. Whether to include an 
#'  ordered factor of plotmath expressions in \code{pg_id_map}
#' @import data.table
#' @return List of two lists (one epithelial; one stromal_immune) 
#'  that contain phenotype characteristics of the 
#'  form 'pg_id_map' - phenotypes table; 'description' - named vector of 
#'  meta_id with description and 'colours' - as for description but
#'  named by colours 
get_phenotypes <- function(TeX_to_expression = F){

  ep_phenos <- here('DATA/Ep_metaclusters.csv')
  str_phenos <- here('DATA/SI_metaclusters.csv')

  get_lab_vec <- function(file){
    
    phenos <- fread(file)
    phenos <- data.table(phenos)

    # Superscripting by Latex syntax for use with latex2exp in plots
    phenos[, meta_description := gsub('MHC\\-', 'MHC ', meta_description)]
    phenos[, meta_description := gsub('med', 'med', meta_description)]
    phenos[, meta_description := gsub('-med', 'med', meta_description)]
    phenos[, meta_description := gsub('-high|-hi|high', 'hi', meta_description)]
    phenos[, meta_description := gsub('-low|-lo|low', 'lo', meta_description)]
    phenos[, meta_description := gsub('(\\+|lo$|med$|hi$|lo |med |hi |-)', '^\\{\\1\\}', meta_description)]
    phenos[, meta_description := gsub(' \\}', '\\}', meta_description)]
    phenos[, meta_description := gsub('reg', 'Reg', meta_description)]
    phenos[, meta_description := gsub('ex', 'Ex', meta_description)]
    phenos[, meta_description := gsub('(Reg|Ex)', '_\\{\\1\\}', meta_description)]
    phenos[, meta_description := gsub('CK8\\^\\{-\\}', 'CK8-', meta_description)]  

    phenos[, meta_id := print_order]
    pg_meta_id_map <- phenos[order(print_order), 
      .(pg_membership, meta_id, meta_description, colours, print_order)]
    meta_description <- unique(pg_meta_id_map$print_order)
    names(meta_description) <- unique(pg_meta_id_map$meta_description)
    colours <- unique(pg_meta_id_map$print_order)
    names(colours) <- unique(pg_meta_id_map$colours)

    if(TeX_to_expression){
      pg_meta_id_map[, meta_expression := factor(print_order, levels = print_order, 
        labels = latex2exp::TeX(meta_description), ordered = TRUE)]
    }

    return(list(
      pg_id_map = pg_meta_id_map,
      description = meta_description,
      colours = colours
      )
    )
  }

  types <- c(ep_phenos, str_phenos)
  out <- lapply(types, get_lab_vec)
  names(out) <- c('epithelial', 'stromal_immune')
  
  cell_labels_correction <-  c('$CD8^{+}$ T cells$','$CD4^{+}$ T cells$','$CD4^{+}$ T cells & APCs$', '$T_{Reg}$ & $T_{Ex}$', 'B cells',
                               '$CD38^{+}$ lymphocytes$','$CD57^{+}$', 'Macrophages', 'Macrophages & granulocytes', 'Granulocytes',
                               '$Fibroblasts$ $FSP1^{+}$','Fibroblasts', 'Myofibroblasts', 'Myofibroblasts $PDPN^{+}$',
                               'Endothelial','$Ki67^{+}$')
  names(out$stromal_immune$description) <- cell_labels_correction
  
  ep_cell_labels_correction <- c("$CK8-18^{hi}CXCL12^{hi}$", "$CK8-18^{+}$ $ER^{hi}$", "$CK8-18^{hi}ER^{lo}$", "$CK^{+}$ $CXCL12^{+}$",
                                 "$ER^{hi}CXCL12^{+}$", "$CK^{med}ER^{lo}$", "$CK^{lo}ER^{med}$", "$CK^{lo}ER^{lo}$", "$HER2^{+}$",
                                 "Basal", "$CD15^{+}$", "$Ep$ $CD57^{+}$", "$MHC^{hi}CD15^{+}$", "$MHC$ $I^{hi}CD57^{+}$",
                                 "$MHC$ $I$ & $II^{hi}$", "$Ep$ K$i67^{+}$")
  names(out$epithelial$description) <- ep_cell_labels_correction
  
  return(out)
}


#' Capitalise first letter; similar to Stata's -proper-
#' 
#' @param string 
#' @return string in proper format
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
      sep="", collapse=" ")
}


#' shorthand for summary
su <- function(x, ...){summary(x, ...)}


#' shorthand for as.data.table
adt <- function(x, ...){as.data.table(x, ...)}


#' shorthand for sort unique elements
sortunique <- function(x, ...){sort(unique(x,...))}


#' Extract substring from character
#'
#'@param pattern
#'@param toMatch character vector
#'@return substring identified by \code[pattern]
getMatched <- function(pattern, toMatch) regmatches(toMatch, regexpr(pattern, toMatch))


#'Load formatted cell neighbour data 
#'
#'@return data.table of formatted cell neighbour data
getCellNeighbours <- function() {
  cells <- read_fst(here('DATA/CellNeighbours.fst'), as.data.table = T)
  return(cells)
}


#'Generate and load spatial epithelial and TME graph data 
#'
#'@return list of two containing three data.tables with graphs and network properties
getSpatialGraphs <- function() {
  source(here::here('Code/Figure4/mkSpatialGraphs.R'), local = TRUE)
  write_fst(communities, here('scratch/communities.fst'))
  saveRDS(NetworkProperties, here('scratch/networkProperties.RDS'))
  return(list(communities = communities, NetworkProperties = NetworkProperties))
}


#'Generate and load patient (tumour) level TME structure data
#'
#'@return data.table of TME structure connectivity data aggregated per structure per tumour
getPtLevelTMEStructures <- function(){
  cellNeighbours <- getCellNeighbours()
  ptLevel <- cellNeighbours[!is.na(TMEStructure), .(ImageNumber, TMEStructure, metabric_id)]
  ptLevel[, TotalInteractions := .N, by = metabric_id]
  ptLevel[,  nInteractionsPerStructure := .N, by = .(metabric_id, TMEStructure)]
  ptLevel <- 
    ptLevel[, .SD[1], by = .(metabric_id, TMEStructure)][, 
    .SD, .SDcols = grep('from|to', names(ptLevel), invert = TRUE, value = TRUE)]
  allCombinations <- adt(expand.grid(ptLevel[, unique(metabric_id)], ptLevel[, unique(TMEStructure)])) #zeros
  setnames(allCombinations, c('Var1', 'Var2'), c('metabric_id', 'TMEStructure'))
  ptLevel <- merge(x = allCombinations, y = ptLevel, by = c('metabric_id', 'TMEStructure'), all.x = TRUE) 
  ptLevel <- ptLevel[order(metabric_id, ImageNumber, decreasing = FALSE)]
  ptLevel[, ImageNumber := ImageNumber[1], by = metabric_id]
  ptLevel[order(TotalInteractions, decreasing = FALSE), TotalInteractions := TotalInteractions[1], by = metabric_id]
  ptLevel[is.na(nInteractionsPerStructure), nInteractionsPerStructure := 0L]
  return(ptLevel)
}


#'Utility returns revised TME structures
#'@return list of data.tables and vectors of TME structure phenotypes
#'
getStructures <- function(){
  file <- here('DATA/TMEStructuresAnnotation.csv')
  dat <- fread(file)
  
  structure_labels_correction <- c("Suppressed expansion", "APC enriched", "TLS-like", "Active IR", "$FSP1^{+}$ enriched$",
                                   "Active stroma", "$PDPN^{+}$ active stroma", "$CD8^{+}$ & $macrophages$", "Granulocyte enriched",
                                   "Vascular stroma")
  dat$Label  <- structure_labels_correction
  
  setkey(dat, 'print_order')
  colours <- dat[, colours]
  names(colours) <- dat[, print_order]
  labs <- dat[, Label]
  names(labs) <- dat[, print_order]


  return(
    list(
      TME = list(
        dat = dat,
        labels = labs,
        colours = colours
      ) 
    ) 
  )
}


#' Get conventional hierarchical clustering order
#' 
#' Uses \code{hclust} to return order 
#' @param X data (df or dt) to cluster
#' @param distance dissim distance
#' @param clustering method for clustering (from option in \code{hclust})
#' @import seriation
#' @return \code{hclust} object (a list that contains an 'order' vector)
getHC <- function(
  X, 
  distance = "euclidean", 
  clustering = "ward.D" 
  )
{

  d <- dist(X, method = distance)
  hc <- hclust(d, method = clustering)
  return(hc)

}



#' FlowSOM warpper for SOM and metaclusters
#' 
#' Makes SOM and, optionally, metaclusters the nodes
#' @param dat dt or df with numeric columns (channels)
#' @param xdim x-dimensions of SOM
#' @param ydim y-dimensions of SOM
#' @param metacluster boolean; should som nodes be metaclustered?
#' @param k  numeric scalar; k of metaclusters if \code{metacluster = T}
#' @param seed integer; random seed to pass to metaClustering_consensus fun
#' @import data.table
#' @return list of som_nodes and, optionally, fsom_clusters

# Define function
getSOM <- function(
  dat, 
  xdim = 10, 
  ydim = 10, 
  metacluster = F, 
  k = 30 , 
  seed = 357951)
{
  
  data_FlowSOM <- flowCore::flowFrame(as.matrix(dat))
  out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
  out <- FlowSOM::BuildSOM(out, xdim = xdim, ydim = ydim)
  som_nodes <- out$map$mapping[, 1]
  
  if(metacluster == T){
    out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = seed)
    fsom_cluster <- out[som_nodes]
    
    FlowSOM_out <- list(
      som_nodes = som_nodes,
      fsom_clusters = fsom_cluster
    )
    
  } else{
    FlowSOM_out <- list(
      som_nodes = som_nodes
    )
  }
  return(FlowSOM_out)
  
}


#' Rescale vector to 0 1
#' 
#' Rescale numeric vector to lie between zero and one, removing \code{NAs}  
#' @param x numeric vector
#' @return numeric vector
rescale01 <- function(x){ 
  x <- (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
  return(x)
}

