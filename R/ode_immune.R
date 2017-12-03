#' Within-host Immune ODE Model
#'
#' write me
#'
#' @param time times to evaluate state variables at
#' @param theta named vector of parameters
#' @param state initial states
#'
#' @export
immune_ode <- function(time, theta, state, ...){
  return(deSolve::ode(y=state,parms=theta,time=time,func="derivs_immune",dllname="RACD",initfunc="init_immune",...))
}
