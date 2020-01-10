# -------------------------------------------------------------------------------- #
#
#   P. falciparum simulation (based on Imperial Model)
#   Equilibrium solutions for mosquito model
#   January 2020
#   Sean Wu (slwu89@berkeley.edu)
#
# -------------------------------------------------------------------------------- #


library(numDeriv)

#' Approximate (continuous-time) Equilibrium for Mosquito Model
#'
#' Based on solutions of the mosquito model as a system of ODEs, return the values of EL, LL, K that
#' hold at equilibrium.
#'
#' @param theta named vector of parameters (see \code{\link{RACD_Parameters}})
#' @param lambdaV equilibrium force of infection on mosquitos
#' @param IV number of infectious vectors at equilibrum
#'
#' @export
RACD_approx_equilibrium <- function(theta,lambdaV,IV){

  x <- rep(0,3)
  x <- setNames(x,c("EL","LL","K"))

  with(as.list(theta),{

    f0 <- tau1+tau2
    a0 <- Q0*f0 # Human biting rate at equilibrium

    EV_eq <- durEV*IV*muV
    SV_eq <- ((durEV*IV*muV) + ((durEV^2)*IV*(muV^2))) / (durEV*lambdaV)
    NV_eq <- IV + EV_eq + SV_eq

    omega_1 <- -0.5 * ( (gamma*(muLL/muEL)) - (durEL/durLL) + ((gamma-1)*muLL*durEL) )
    omega_2 <- sqrt(
      (0.25 * (( (gamma*(muLL/muEL)) - (durEL/durLL) + ((gamma-1)*muLL*durEL) )^2)) + (gamma * ((beta*muLL*durEL) / (2*muEL*muV*durLL*(1 + (durPL*muPL)))) )
    )
    omega <- omega_1 + omega_2

    PL_eq <- 2*durPL*muV*NV_eq
    LL_eq <- 2*muV*durLL*(1 + (durPL*muPL))*NV_eq
    EL_eq <- 2*omega*muV*durLL*(1 + (durPL*muPL))*NV_eq

    K <- (NV_eq*2*durLL*muV*(1 + (durPL*muPL))*gamma*(omega+1)) / ((omega/(muLL*durEL)) - (1/(muLL*durLL)) - 1)

    x[["EL"]] <<- EL_eq
    x[["LL"]] <<- LL_eq
    x[["K"]] <<- K
  })

  return(x)
}

#' Objective Function for Mosquito Equilibrium
#'
#' Calculate finite differences of discrete-time model
#'
#' @param theta named vector of parameters (see \code{\link{RACD_Parameters}})
#' @param dt size of time-step
#' @param NV total size of vector population
#' @param PL size of pupae population
#'
#'
make_obj_f <- function(theta,dt,NV,PL){

  # parameters
  NV <- NV
  dt <- dt
  beta <- theta[["beta"]]
  muEL <- theta[["muEL"]]
  durEL <- theta[["durEL"]]
  muLL <- theta[["muLL"]]
  gamma <- theta[["gamma"]]
  durLL <- theta[["durLL"]]
  muPL <- theta[["muPL"]]
  durPL <- theta[["durPL"]]
  PL <- PL

  fn <- function(x){

    # variables to optimize
    EL <- x[1]
    LL <- x[2]
    K <- x[3]

    EL_mort <- muEL*(1 + ((EL+LL)/K))
    r_EL <- 1/durEL

    deltaE <- (beta*NV*dt) - ((1 - exp(-(EL_mort + r_EL)*dt)) * EL)

    LL_mort <-muLL*(1 + gamma*((EL+LL)/K))
    r_LL <- 1/durLL

    deltaL <- ((1 - exp(-(EL_mort + r_EL)*dt)) * (r_EL / (EL_mort + r_EL)) * EL) - ((1 - exp(-(LL_mort + r_LL)*dt)) * LL)

    r_PL <- 1/durPL

    deltaP <- ((1 - exp(-(LL_mort + r_LL)*dt)) * (r_LL / (LL_mort + r_LL)) * LL) - ((1 - exp(-(muPL + r_PL)*dt)) * PL)

    return(sum(deltaE^2,deltaL^2,deltaP^2))
  }

  return(fn)
}


#' Linear Constraints Function for Mosquito Equilibrium
#'
#' Calculate linear constraints
#'
#' @param x numeric vector
#'
cons_g <- function(x){
  EL <- x[1]
  LL <- x[2]
  K <- x[3]
  rbind(
    EL,
    LL,
    K,
    EL - LL
  )
}

#' Numeric Jacobian of Linear Constraints Function for Mosquito Equilibrium
#'
#' Calculate Jacobian matrix (1st order partial derivatives) of linear constraints
#'
#' @param x numeric vector
#'
jac_g <- function(x){
  numDeriv::jacobian(cons_g,x)
}


#' Exact (discrete-time) Equilibrium for Mosquito Model
#'
#' Using approximate solutions from \code{\link[RACD]{RACD_approx_equilibrium}} as starting point(s)
#' use constrained numerical optimization to find equilibrium solutions of (EL, LL, K)
#' from a given force of infection on mosquitoes and size of infected vector population,
#' which are obtainable from solving the human model at equilibrium.
#' Starting values for the optimzation are found by using Latin Hypercube Sampling
#' in a box around the approximate solutions.
#'
#' @param theta named vector of parameters (see \code{\link{RACD_Parameters}})
#' @param dt size of time-step
#' @param IV number of infectious vectors at equilibrum
#' @param lambdaV equilibrium force of infection on mosquitos
#' @param sd size of box around approximate value is (x +- x*sd)
#' @param nstart number of starting points for optimization
#'
#' @export
RACD_mosq_equilibrium <- function(theta,dt,IV,lambdaV,sd=0.25,nstart=250){

  # get approximate values from the continuous-time model
  approx <- RACD_approx_equilibrium(theta,lambdaV,IV)

  # (some of) the equilibrium values for state variables
  SV <- EV <- PL <- NV <- 0

  with(as.list(theta),{

    # jump probabilities
    IV2D <- 1 - exp(-muV*dt)
    EV2D <- (1 - exp(-((1/durEV) + muV)*dt)) * (muV / ((1/durEV) + muV))
    EV2IV <- (1 - exp(-((1/durEV) + muV)*dt)) * ((1/durEV) / ((1/durEV) + muV))
    SV2EV <- (1 - exp(-(lambdaV + muV)*dt)) * (lambdaV / (lambdaV + muV))
    # equilibria for SV and EV
    SV <<- (IV*IV2D + ((EV2D*IV*IV2D) / EV2IV)) / SV2EV
    EV <<- (IV*IV2D)/EV2IV

    # jump probabilities
    SV2ALL <- (1 - exp(-(lambdaV + muV)*dt))
    P2SV <- (1 - exp(-(muPL + (1/durPL))*dt)) * ((1/durPL) / ((1/durPL) + muPL))
    # equilibria for P
    PL <<- (2*SV*SV2ALL) / P2SV

    # equilibria total vector population
    NV <<- SV + EV + IV

  })

  # use numerical optimization to handle the remaining highly non-linear equations
  # (several hours of my own time and mathematica did not produce anything reasonable analytically)
  # inequality constraints:
  # EL > LL
  # positivity for all variables

  # objective fn: minimize SSE
  obj_f <- make_obj_f(theta,dt,NV,PL)

  grad_f <- function(x){
    numDeriv::grad(obj_f, x, method="Richardson")
  }

  if(any(cons_g(approx) < 0)){
    stop("initial values do not respect constraints!")
  }

  # linear constraints for constrOptim (k x p)
  ui <- matrix(c(1,0,0,
                 0,1,0,
                 0,0,1,
                 1,-1,0),
               byrow = T,nrow = 4,ncol = 3)
  ci <- rep(0,4) # k

  # from many starting points
  mins <- approx - (approx*sd)
  maxs <- approx + (approx*sd)
  p <- ncol(ui)
  starts <- as.list(data.frame(t(lhs::randomLHS(n = nstart,k = p))))
  starts <- lapply(starts,function(x){
    qunif(p = x,min = mins,max = maxs)
  })

  control <- list(trace=0,maxit=1e3,reltol=1e-10)
  opt <- lapply(X = starts,FUN = function(start){
    constrOptim(theta = start,f = obj_f,grad = grad_f,
                ui = ui,ci = ci, method = "Nelder-Mead",
                control=control,outer.eps = 1e-6,outer.iterations = 2e2)
  })

  min <- which.min(sapply(opt,function(x){x$value}))

  # return analytic and numerically approximated solutions
  out <- list(
    SV_eq = SV,
    EV_eq = EV,
    IV_eq = IV,
    PL_eq = PL,
    LL_eq = opt[[min]]$par[2],
    EL_eq = opt[[min]]$par[1],
    K_eq = opt[[min]]$par[3],
    NV_eq = NV,
    opt_min = opt[[min]]
  )
  return(out)
}
