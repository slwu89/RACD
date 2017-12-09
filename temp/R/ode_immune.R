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

  immune = RACD::immune_ode(time = time_immune,theta = theta_immune,state = state_immune,...)
  return(c(IB=unname(immune[nrow(immune),"IB"]),ID=unname(immune[nrow(immune),"ID"]),ICA=unname(immune[nrow(immune),"ICA"])))
}


###############################################################################
# Within-host Infection ODE Model
###############################################################################

#' Within-host Infection ODE Model
#'
#' write me
#'
#' @param time times to evaluate state variables at
#' @param theta named vector of parameters
#' @param state initial states
#' @param ... additional parameters passed to numerical integrator \code{\link[deSolve]{ode}}
#'
#' @export
infection_ode <- function(time, theta, state, ...){
  return(deSolve::ode(y=state,parms=theta,time=time,func="derivs_infection",dllname="RACD",initfunc="init_infection",...))
}

#' Within-host Infection State ODE Model
#'
#' ODE model for calculating the probability that an individual is in each state given their age and EIR heterogeneity attributes.
#' This function calls \code{\link{infection_ode}}
#'
#' @param a age of individual
#' @param zeta individual biting heterogeneity
#' @param psi relative biting rate (used to calculate psi(a))
#' @param theta vector of model parameters (see \code{\link{RACD_Parameters}})
#' @param ... additional parameters passed to numerical integrator \code{\link[deSolve]{ode}}
#'
#' @export
infection_initial <- function(a, zeta, psi, theta, ...){

  time_infection = seq(from=0,to=a,by=a/1e3)

  init_immune = RACD::immune_initial(a=a,zeta=zeta,psi=psi,theta=theta,...)

  initICA20 = init_immune[["ICA"]]

  state_infection = c(prS=1, prT=0, prD=0, prA=0, prU=0, prP=0,IB=0, ICA=0)

  theta_infection = c(theta[c("epsilon0","fT","dE","dT","dD","dA","dU","dP","rho","a0","b0","b1","dB","IB0","kappaB","uB","phi0","phi1","dC","IC0","kappaC","uC","PM","dM")],"initICA20"=initICA20,"zeta"=zeta,"psi"=psi)
  theta_infection[["epsilon0"]] = theta_infection[["epsilon0"]] * 365
  theta_infection[["dE"]] = theta_infection[["dE"]] / 365
  theta_infection[["dT"]] = theta_infection[["dT"]] / 365
  theta_infection[["dD"]] = theta_infection[["dD"]] / 365
  theta_infection[["dA"]] = theta_infection[["dA"]] / 365
  theta_infection[["dU"]] = theta_infection[["dU"]] / 365
  theta_infection[["dP"]] = theta_infection[["dP"]] / 365
  theta_infection[["dB"]] = theta_infection[["dB"]] / 365
  theta_infection[["uB"]] = theta_infection[["uB"]] / 365
  theta_infection[["dC"]] = theta_infection[["dC"]] / 365
  theta_infection[["uC"]] = theta_infection[["uC"]] / 365
  theta_infection[["dM"]] = theta_infection[["dM"]] / 365

  infection = RACD::infection_ode(time_infection,theta_infection,state_infection,...)
  return(c(
    prS = infection[nrow(infection),"prS"],
    prT = infection[nrow(infection),"prT"],
    prD = infection[nrow(infection),"prD"],
    prA = infection[nrow(infection),"prA"],
    prU = infection[nrow(infection),"prU"],
    prP = infection[nrow(infection),"prP"]
  ))
}
