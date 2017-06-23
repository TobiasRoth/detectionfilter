#' Two surveys of Swiss plant communities from Swiss Biodiversity Monitoring.
#' 
#' The Swiss biodiversity monitoring (BDM) has among others monitored the 
#' presence/absence of plant species at 1x1km squares. The data-set contains the
#' data from 363 squares with exactly two surveys conducted by professional 
#' botanists between 2010 and 2014.
#' 
#' 
#' @format A list with xx elements: \describe{ \item{y}{Array with 
#'   presence/absence observed in 362 squares for 1733 species during two 
#'   surveys.} 
#'   \item{elevation}{Vecotor with elevation (m) of each square.} 
#'   \item{dates}{Data.frame containing the Julian days for the two visits that
#'   were conducted to each square. The line correspond to the squares and the
#'   two columns to the two surveys.}
#'   
#'   
#'   }
#' @source \url{http://www.biodiversitymonitoring.ch/}
"plantsBDM"


#' Traits for the 1733 species observed during the BDM surveys.
#' 
#' The Swiss 
#' 
#' 
#' @format A data.frame with the following collumns: 
#' \describe{ 
#'   \item{sla}{Specific leaf area.} 
#'   \item{ch}{Canopy height.} 
#'   \item{sm}{Seed mass.}
#'   
#'   
#'   }
#' @source \url{http://www.uni-oldenburg.de/en/biology/landeco/research/projects/leda/}
"traitmat"




