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

#' Process RACD Simulation Output
#'
#' Given output from \code{\link[RACD]{RACD_Simulation}}, take discrete event output
#' and make daily state occupancy vector as matrix.
#'
#' @param outfile path to logged events
#'
#' @export
RACD_StateVector <- function(outfile){

  # read in raw data
  raw = read.table(file = outfile,header = TRUE,sep = ",",stringsAsFactors = FALSE)

  # set up state occupancy vector
  minT = min(raw$Time)
  maxT = max(raw$Time)
  states = c("S","E","T","D","A","U","P")
  state_vector = matrix(0L,nrow=(maxT-minT)+1,ncol=length(states)+1,dimnames = list(NULL,c("time",states)))
  state_vector[,"time"] = minT:maxT

  # humans
  humans = unique(raw$HumanID)
  pb = txtProgressBar(min = 0,max = max(humans),style = 3)
  for(h in humans){
    hh = raw[raw$HumanID==h,]
    
    # single event
    if(nrow(hh)==1){
      if(hh$Event=="Birth"){
        state_vector[hh$Time:nrow(state_vector),"S"] = state_vector[hh$Time:nrow(state_vector),"S"]+ 1L
      } else {
        state_vector[hh$Time:nrow(state_vector),hh$Event] = state_vector[hh$Time:nrow(state_vector),hh$Event] + 1L
      }
    } else {
      # multiple events: handle up to last event
      for(e in 1:(nrow(hh)-1)){
        if(hh$Event[e]=="Birth"){
          state_vector[hh$Time[e]:(hh$Time[e+1]-1),"S"] = state_vector[hh$Time[e]:(hh$Time[e+1]-1),"S"] + 1L
        } else {
          state_vector[hh$Time[e]:(hh$Time[e+1]-1),hh$Event[e]] = state_vector[hh$Time[e]:(hh$Time[e+1]-1),hh$Event[e]] + 1L
        }
      }
      # multiple events: handle last event if not death
      if(hh$Event[e+1]!="Death"){
        state_vector[hh$Time[e+1]:nrow(state_vector),hh$Event[e+1]] = state_vector[hh$Time[e+1]:nrow(state_vector),hh$Event[e+1]] + 1L 
      }
    }
    
    setTxtProgressBar(pb,value = h)
  }
  
  return(state_vector)
}
