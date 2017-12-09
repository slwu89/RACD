###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   December 2017
#   ODE Immune Model Wrappers
#
###############################################################################


###############################################################################
# Within-host Immune ODE Model
###############################################################################

#' Within-host Immune ODE Model
#'
#' write me
#'
#' @param time times to evaluate state variables at
#' @param theta named vector of parameters
#' @param state initial states
#' @param ... additional parameters passed to numerical integrator \code{\link[deSolve]{ode}}
#'
#' @export
immune_ode <- function(time, theta, state, ...){
  return(deSolve::ode(y=state,parms=theta,time=time,func="derivs_immune",dllname="RACDaux",initfunc="init_immune",...))
}

#' Within-host Pre-erythrocytic Immune ODE Model
#'
#' Calculate an individual's pre-erythrocytic immunity at age a for given EIR heterogeneities (biting rate & geographic location).
#' This function calls \code{\link{immune_ode}}
#'
#' @param a age of individual
#' @param zeta individual biting heterogeneity
#' @param psi relative biting rate (used to calculate psi(a))
#' @param theta vector of model parameters (see \code{\link{RACD_Parameters}})
#' @param ... additional parameters passed to numerical integrator \code{\link[deSolve]{ode}}
#'
#' @export
immune_initial <- function(a, zeta, psi, theta, ...){

  theta_immune = c("zeta"=zeta,"psi"=psi,theta[c("a0","rho","dB","uB","dID","uD","dC","uC","b0","b1","IB0","kappaB","epsilon0")])
  theta_immune[["dB"]] = theta_immune[["dB"]]/365
  theta_immune[["uB"]] = theta_immune[["uB"]]/365
  theta_immune[["dID"]] = theta_immune[["dID"]]/365
  theta_immune[["uD"]] = theta_immune[["uD"]]/365
  theta_immune[["dC"]] = theta_immune[["dC"]]/365
  theta_immune[["uC"]] = theta_immune[["uC"]]/365
  theta_immune[["epsilon0"]] = theta_immune[["epsilon0"]]*365

  state_immune = c("IB"=0,"ID"=0,"ICA"=0)

  time_immune = seq(from=0,to=a,by=a/1e3)

  immune = RACDaux::immune_ode(time = time_immune,theta = theta_immune,state = state_immune,...)
  return(c(
    IB = unname(immune[nrow(immune),"IB"]),
    ID = unname(immune[nrow(immune),"ID"]),
    ICA = unname(immune[nrow(immune),"ICA"])
  ))
}
