###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   December 2017
#   Process Output
#
###############################################################################

#' #' Process RACD Simulation Output
#' #'
#' #' Given output from \code{\link[RACD]{RACD_Simulation}}, take continuous time state transition output
#' #' and make state occupancy vector by binning events by day.
#' #'
#' #' @param outfile path to logged events
#' #'
#' #' @export
#' RACD_StateVector <- function(outfile){
#' 
#'   # read in raw data
#'   raw = read.table(file = outfile,header = TRUE,sep = ",",stringsAsFactors = FALSE)
#' 
#'   # set up state occupancy vector
#'   minT = min(raw$Time)
#'   maxT = max(raw$Time)
#'   states = c("S","E","T","D","A","U","P")
#'   state_vector = matrix(0,nrow=(maxT-minT)+1,ncol=length(states)+1,dimnames = list(NULL,c("time",states)))
#'   state_vector[,"time"] = minT:maxT
#' 
#'   # calculate initial distribution of states
#'   # t = state_vector[i,"time"]
#'   # tab = table(raw[raw$Time==t,"Event"])
#'   # state_vector[i,match(names(tab),colnames(state_vector))] = tab
#' 
#'   # humans
#'   humans = unique(raw$HumanID)
#'   for(h in humans){
#'     hh = raw[raw$HumanID==h,]
#'   }
#' 
#' }
