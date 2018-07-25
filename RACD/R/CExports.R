# This was written by hand
# This is probably not the proper way to do this
# Do not use this as an example

#' @useDynLib RACD init_immune
something0 <- function(integer0, double0) .Call(init_immune, integer0, double0)

#' @useDynLib RACD derivs_immune
something1 <- function(integer0, double0, double1, double2, doubl3, integer1){
  .Call(derivs_immune, integer0, double0, double1, double2, doubl3, integer1)
}

#' @useDynLib RACD init_infection
something2 <- function(integer0, double0) .Call(init_infection, integer0, double0)
  
#' @useDynLib RACD derivs_infection
something3 <- function(integer0, double0, double1, double2, doubl3, integer1){
  .Call(derivs_infection, integer0, double0, double1, double2, doubl3, integer1)
}
