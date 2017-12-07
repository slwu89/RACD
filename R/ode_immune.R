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
#' @param ... additional parameters passed to numerical integrator \code{\link[deSolve]{ode}}
#'
#' @export
immune_initial <- function(a, zeta, psi,
  a0=8,rho=0.85,
  durB=(3650/365),uB=(7.2/365),b0=0.59,b1=0.5,IB0=43.9,kappaB=2.16,
  durD=(3650/365),uD=(9.45/365),
  durC=(10950/365),uC=(6.06/365),
  epsilon0=(0.01369863*365), ...
){
  times = seq(from=0,to=a,by=a/1000)
  theta =  c(a0=a0,rho=rho,zeta=zeta,psi=psi,durB=durB,uB=uB,b0=b0,b1=b1,IB0=IB0,kappaB=kappaB,durD=durD,uD=uD,durC=durC,uC=uC,epsilon0=epsilon0)
  immune = RACD::immune_ode(time = times,theta = theta,state = c(IB=0, ID=0, ICA=0),...)
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
#'
#' @param a age of individual
#' @param zeta individual biting heterogeneity
#' @param psi relative biting rate (used to calculate psi(a))
#' @param epsilon0 mean EIR
#' @param fT proportion of clinical disease cases successfully treated
#' @param dE Duration of latent period (years)
#' @param dT Duration of treated clinical disease (years)
#' @param dD Duration of untreated clinical disease (years)
#' @param dA Duration of patent infection (years)
#' @param dU Duration of sub-patent infection (years) (fitted)
#' @param dP Duration of prophylactic protection following treatment (years)
#' @param rho age-dependent biting heterogeneity (used to calculate psi(a))
#' @param a0 age-dependent biting heterogeneity (used to calculate psi(a))
#' @param b0 probability of infection with no immunity
#' @param b1 maximum relative reduction in immunity
#' @param dB Inverse of decay rate (years) of immunity reducing probability of infection
#' @param IB0 scale parameter for probability of infection
#' @param kappaB shape parameter for probability of infection
#' @param durB durB 1/decay rate of immunity reducing probability of infection
#' @param uB duration in which immunity is not boosted
#' @param durD 1/decay rate of immunity reducing probability of detection
#' @param uD duration in which immunity is not boosted
#' @param durC 1/decay rate of immunity reducing probabiltiy of clinical disease
#' @param uC duration in which immunity is not boosted
#' @param phi0 probability of clinical disease with no immunity
#' @param phi1 maximum relative reduction of clinical disease for immunity reducing probability of clinical disease
#' @param dC 1/decay rate of immunity reducing probability of clinical disease
#' @param IC0 Scale parameter of immunity reducing probability of clinical disease
#' @param kappaC Shape parameter of immunity reducing probability of clinical disease
#' @param PM New-born immunity relative to mother's immunity
#' @param dM Inverse decay rate of maternal immunity (years)
#' @param ... additional parameters passed to numerical integrator \code{\link[deSolve]{ode}}
#'
#' @export
infection_initial <- function(a, zeta, psi,
  epsilon0 = (0.01369863*365), fT = 0.4,
  dE = (12/365), dT = (5/365), dD = (5/365), dA = (200/365), dU = (110/365), dP = (25/365),
  rho = 0.85, a0 = 8,
  b0 = 0.59, b1 = 0.5, dB = (3650/365), IB0 = 43.9, kappaB = 2.16,
  durB=(3650/365), uB = (7.2/365),
  durD=(3650/365),uD=(9.45/365),
  durC=(10950/365),uC=(6.06/365),
  phi0 = 0.792, phi1 = 0.00074, dC = (10950/365), IC0 = 18, kappaC = 2.37, PM = 0.774, dM = (67.7/365), ...
){
  initI20 = immune_initial(a=a, zeta=zeta, psi=psi, a0=a0,rho=rho,durB=durB,uB=uB,b0=b0,b1=b1,IB0=IB0,kappaB=kappaB,durD=durD,uD=uD,durC=durC,uC=uC,epsilon0=epsilon0,...)

  initICA20 = initI20[["ICA"]]
  times = seq(from=0,to=a,by=a/1000)
  stateInf = c(prS=1, prT=0, prD=0, prA=0, prU=0, prP=0,IB=0, ICA=0)
  thetaInf = c(epsilon0=epsilon0,fT=fT,dE=dE,dT=dT,dD=dD,dA=dA,dU=dU,dP=dP,rho=rho,a0=a0,b0=b0,b1=b1,dB=dB,IB0=IB0,kappaB=kappaB,uB=uB,phi0=phi0,phi1=phi1,dC=dC,IC0=IC0,kappaC=kappaC,uC=uC,PM=PM,dM=dM,initICA20=initICA20,zeta=zeta,psi=psi)
  initImm = RACD::infection_ode(time=times,theta=thetaInf,state=stateInf,...)
  return(c(
    prS = initImm[nrow(initImm),"prS"],
    prT = initImm[nrow(initImm),"prT"],
    prD = initImm[nrow(initImm),"prD"],
    prA = initImm[nrow(initImm),"prA"],
    prU = initImm[nrow(initImm),"prU"],
    prP = initImm[nrow(initImm),"prP"]
  ))
}
