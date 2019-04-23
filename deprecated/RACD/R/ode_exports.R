#' @useDynLib RACD init_immune
export_init_immune <- function(integer0, double0) .Call(init_immune, integer0, double0)

#' @useDynLib RACD derivs_immune
export_derivs_immune <- function(integer0, double0, double1, double2, doubl3, integer1){
  .Call(derivs_immune, integer0, double0, double1, double2, doubl3, integer1)
}

#' @useDynLib RACD init_infection
export_init_infection <- function(integer0, double0) .Call(init_infection, integer0, double0)

#' @useDynLib RACD derivs_infection
export_derivs_infection <- function(integer0, double0, double1, double2, doubl3, integer1){
  .Call(derivs_infection, integer0, double0, double1, double2, doubl3, integer1)
}
