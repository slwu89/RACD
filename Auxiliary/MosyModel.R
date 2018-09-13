###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   September 2018
#   Tests of Mosquito Model
#
###############################################################################


# PHD system for adults

# make Q (rate matrix)
#' @param
makeQ <- function(lambda,q,n,mu){

  # rate matrix
  D <- 3+n
  Q <- Matrix::Matrix(data=0,nrow=D,ncol=D)

  # S vectors
  Q[1,2] <- lambda
  Q[1,D] <- mu
  Q[1,1] <- -sum(Q[1,])

  # E vectors
  Q[2,3] <- q*n
  Q[2,D] <- mu
  Q[2,2] <- -sum(Q[2,])

  for(i in (2:n)+1){
    Q[i,i+1] <- q*n
    Q[i,D] <- mu
    Q[i,i] <- -sum(Q[i,])
  }
  i <- i + 1

  # I vectors
  Q[i,D] <- mu
  Q[i,i] <- -sum(Q[i,])

  dimnames(Q) <- list(
    c("S",paste0("E",1:n),"I","D"),
    c("S",paste0("E",1:n),"I","D")
  )

  return(Q)
}

makeQ_char <- function(lambda,q,n,mu){

  # rate matrix
  D <- 3+n
  Q <- matrix(data="0",nrow=D,ncol=D)

  # S vectors
  Q[1,2] <- "lambda"
  Q[1,D] <- "mu"
  Q[1,1] <- "-(lambda + mu)"

  # E vectors
  Q[2,3] <- "q*n"
  Q[2,D] <- "mu"
  Q[2,2] <- "-(q*n + mu)"

  for(i in (2:n)+1){
    Q[i,i+1] <- "q*n"
    Q[i,D] <- "mu"
    Q[i,i] <- "-(q*n + mu)"
  }
  i <- i + 1

  # I vectors
  Q[i,D] <- "mu"
  Q[i,i] <- "-mu"

  dimnames(Q) <- list(
    c("S",paste0("E",1:n),"I","D"),
    c("S",paste0("E",1:n),"I","D")
  )

  return(Q)
}

# ODE system
param <- c(
  beta = 21.19, # Number of eggs laid per day by female mosquito
  muEL = 0.034, # Early larval instar daily mortality
  muLL = 0.035, # Late larval instar daily mortality
  muPL = 0.25, # Pupal daily mortality
  durEL = 6.64, # Duration of early instar stage
  durLL = 3.72, # Duration of late instar stage
  durPL = 0.64, # Duration of pupal stage
  gamma = 13.25, # Effect of density-dependence on late instars relative to early instars
  muV = 1/7.6, # Adult mosquito daily mortality
  q = 1/10, # mean length of EIP
  n = 64, # shape parameter for Gamma distributed EIP
  f0 = 1/3, # Daily biting rate by mosquitoes on animals and humans
  Q0 = 0.92 # human blood index
)

epsilon0 = 10/365, # Daily entomological inolculation rate
iH_eq = 0.35, # Equilibrium malaria prevalence in humans
NH_eq = 200, # Equilibrium human population size
cV = 0.05 # Probability of transmission from human to vector per infectious bite

# derived parameters
with(as.list(param),{
  b_omega <<- gamma*muLL/muEL - durEL/durLL + (gamma-1)*muLL*durEL
  omega <<- -0.5*b_omega + sqrt(0.25*b_omega^2 + gamma*beta*muLL*durEL/(2*muEL*muV*durLL*(1+durPL*muPL)))
  a0 <<- Q0*f0 # Human biting rate at equilibrium
})
