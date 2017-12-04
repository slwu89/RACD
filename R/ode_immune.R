###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   December 2017
#   ODE Model Wrappers
#
###############################################################################

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

#' Within-host Pre-erythrocytic Immune ODE Model
#'
#' Calculate an individual's pre-erythrocytic immunity at age a for given EIR heterogeneities (biting rate & geographic location).
#' This function calls \code{\link{immune_ode}}
#'
#' @param a age of individual
#' @param zeta individual biting heterogeneity
#' @param psi relative biting rate (used to calculate psi(a))
#' @param a0 age-dependent biting heterogeneity (used to calculate psi(a))
#' @param rho age-dependent biting heterogeneity (used to calculate psi(a))
#' @param durB 1/decay rate of immunity reducing probability of infection
#' @param uB duration in which immunity is not boosted
#' @param b0 probability of infection with no immunity
#' @param b1 maximum relative reduction in immunity
#' @param IB0 scale parameter for probability of infection
#' @param kappaB shape parameter for probability of infection
#' @param durD 1/decay rate of immunity reducing probability of detection
#' @param uD duration in which immunity is not boosted
#' @param durC 1/decay rate of immunity reducing probabiltiy of clinical disease
#' @param uC duration in which immunity is not boosted
#' @param epsilon0 mean EIR
#'
#' @export
immune_initial <- function(a, zeta, psi,
           a0=8,rho=0.85,
           durB=(3650/365),uB=(7.2/365),b0=0.59,b1=0.5,IB0=43.9,kappaB=2.16,
           durD=(3650/365),uD=(9.45/365),
           durC=(10950/365),uC=(6.06/365),epsilon0=(0.01369863*365),...){

  times = seq(from=0,to=a,by=a/1000)
  theta = c(a0=a0,rho=rho,zeta=zeta,psi=psi,durB=durB,uB=uB,b0=b0,b1=b1,IB0=IB0,kappaB=kappaB,durD=durD,uD=uD,durC=durC,uC=uC,epsilon0=epsilon0)
  immune = RACD::immune_ode(time = times,theta = theta,state = c(IB=0, ID=0, ICA=0),...)
  return(c(IB=immune[nrow(immune),"IB"],ID=immune[nrow(immune),"ID"],ICA=immune[nrow(immune),"ICA"]))
}

#' Within-host Infection ODE Model
#'
#' write me
#'
#' @param time times to evaluate state variables at
#' @param theta named vector of parameters
#' @param state initial states
#'
#' @export
infection_ode <- function(time, theta, state, ...){
  return(deSolve::ode(y=state,parms=theta,time=time,func="derivs_infection",dllname="RACD",initfunc="init_infection",...))
}
