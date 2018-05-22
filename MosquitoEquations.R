###############################################################################
###############################################################################
# Malaria vector ODE model
# Jared Bennett
# 2018/05/22
###############################################################################
###############################################################################
# This is the complete set of life-cycle and equilibrium solutions for the mosquito
#  model in RACD.

###############################################################################
# Initialization
###############################################################################
# Given a shit-ton of parameters, the following equations will solve for the
#  equilibrium solution to initialize a mosquito population.

#############################
# Parameters
#############################
Parameters1 <- c(
  # Human parameters
  "EIR" = ,
  "H" = ,
  "I_H" = ,

  # transmission parameters
  "a" = ,
  "b_V" = ,
  "tau" = ,

  # Life parameters
  "d_E" = ,
  "d_L" = ,
  "d_P" = ,
  "mu_E" = ,
  "mu_P" = ,
  "mu_F" = ,
  "beta" = ,
  "gamma" =
)

#############################
# Equations
#############################
calc_omega <- function(Parameters1){
  # inputs
  #   mu_E, mu_L, mu_P, d_E, d_L, d_P, gamma, beta
  with(Parameters1, {

    term1 <- -1/2*(gamma*mu_L/mu_E - d_E/d_L + (gamma-1)*mu_L*d_E)

    term2 <- 1/4(gamma*mu_L/mu_E - d_E/d_L + (gamma-1)*mu_L*d_E)^2

    term3 <- gamma*beta*mu_L*d_E/(2(mu_E*mu_F*d_L*(1+d_P*mu_P)))

    return(term1 + sqrt(term2 + term3))

  }) # end with
}

calc_MosquitoEquilibrium <- function(Parameters1){
  # inputs
  #   EIR, H, I_H, a, b_V, tau, d_L, d_P, mu_P, mu_F, omega
  with(Parameters1, {

    #E
    numerator <- 2*d_L*mu_F*omega*(1+d_P*mu_P)*(a*b_V*I_H + mu_F*H) * H * exp(mu_F*tau)
    denominator <- a^2*b_V*I_H
    E_eq <- numerator / denominator * EIR

    #L
    numerator <- 2*d_L*mu_F*(1+d_P*mu_P)*(a*b_V*I_H + mu_F*H) * H * exp(mu_F*tau)
    L_eq <- numerator / denominator * EIR

    #P
    numerator <- 2*d_P*mu_F*(a*b_V*I_H + mu_F*H) * H * exp(mu_F*tau)
    L_eq <- numerator / denominator * EIR

    #S_V
    numerator <- mu_F * H^2 * exp(mu_F * tau)
    S_Veq <- numerator / denominator * EIR

    #E_V
    numerator <- mu_F * H * (exp(mu_F * tau)-1)
    E_Veq <- numerator / a * EIR

    #I_V
    I_Veq <- H / a * EIR

    return(c("E" = E_eq,
             "L" = L_eq,
             "P" = P_eq,
             "S_V" = S_Veq,
             "E_V" = E_Veq,
             "I_V" = I_Veq
    )) # end return
  }) # end with
}

calc_K <- function(Parameters1, MosEQ){
  # inputs
  #   MosEQ, d_E, d_L, d_P, mu_P, mu_F, gamma, omega
  with(as.list(c(Parameters1, MosEQ)),{

    F_eq <- sum(MosEQ[c("S_V", "E_V", "I_V")])
    numerator <- 2*d_L*mu_F*(1+d_P*mu_P)*gamma*(omega+1)
    denominator <- omega/(mu_L*d_E) - 1/(mu_L*d_L) - 1

    return(numerator / denominator * F_eq)
  }) # end with
}

###############################################################################
# Interventions
###############################################################################
# Given a just "few" more parameters, and the initialization products from above,
#  the following DDE's describe the life cycle and malaria issues of our mosquitoes.

#############################
# Parameters
#############################
Parameters2 <- c(
  # aquatic interventions
  "nu" = ,
  "rho" = ,

  # adult interventions
  "mu_1" = ,
  "mu_2" = ,
  "t_1" = ,
  "t_2" = ,

  "chi_ITN" = ,
  "chi_IRS" = ,
  "phi_b" = ,

  "r_ITN" = ,
  "r_IRS" = ,

  "sigma_O" = ,
  "sigma_ITN" = ,
  "sigma_IRS" =
)

#############################
# Equations
#############################
calc_Pr <- function(){
  # inputs
  #   chi_ITN, chi_IRS, phi_b, r_ITN, r_IRS

  term1 <- (chi_ITN - chi_ITN*chi_IRS) * phi_b * r_ITN

  term2 <- (chi_IRS - chi_ITN*chi_IRS) * (1 - phi_b) * r_IRS

  term3 <- chi_ITN*chi_IRS * ((1-phi_b)*r_IRS + phi_b*(r_IRS + (1-r_IRS)*r_ITN))

  return(term1 + term2 + term3)
}

calc_Pf <- function(){
  # inputs
  #   chi_ITN, chi_IRS, phi_b, r_ITN, r_IRS, sigma_O, sigma_ITN, sigma_IRS

  term1 <- (1 - chi_ITN - chi_IRS + chi_ITN*chi_IRS) * sigma_O

  term2 <- (chi_ITN - chi_ITN*chi_IRS) * phi_b * (1 - r_ITN) * sigma_ITN

  term3 <- (chi_IRS - chi_ITN*chi_IRS) * (1 - phi_b) * (1 - r_IRS) * sigma_IRS

  term4 <- chi_ITN*chi_IRS * ((1-phi_b)*(1-r_IRS)*sigma_IRS +
                                phi_b*(1 - (r_IRS+(1-r_IRS)*r_ITN)) * sigma_ITN*sigma_IRS)

  return(term1 + term2 + term3 + term4)
}

calc_muFTheta <- function(){
  # inputs
  #   mu_1, mu_2, t_1, t_2, Pr, Pf

  term1_numerator <- mu_1*t_1 + mu_2*t_2*(1-Pr)
  term1_denominator <- t_1 + t_1*(1-Pr)

  term2_base <- log( Pf/(1-Pr*exp( -mu_1*t_1/(1-Pr) )) )
  term2_exp <- (1-Pr)/(t_1 + t_2*(1-Pr))

  return(term1_numerator / term1_denominator - term2_base^term2_exp)
}

###############################################################################
# Life Cycle Equations
###############################################################################
calc_dE <- function(){
  # inputs
  #  beta, S_V, E_V, I_V, E, L, nu, mu_E, K, d_E

  term1_2 <- beta*(S_V + E_V+I_V) - E/d_E

  term3 <- nu*mu_E*(1 + (E+L)/K)*E

  return(term1_2 - term3)
}

calc_dL <- function(){
  # inputs
  #   E, L, K, d_E, d_L, nu, mu_L, gamma

  terms1_2 <- E/d_E - L/d_L

  term3 <- nu*mu_L*(1 + gamma*(E+L)/K)*L

  return(terms1_2 - term3)
}

calc_dP <- function(){
  # inputs
  #   L, P, d_L, nu, mu_P, d_P

  term2 <- nu*(mu_P + 1/d_P)*P

  return(L/d_L - term2)
}

calc_dS_V <- function(){
  # inputs
  #   P, I_H, H, S_V, b_V, a, rho, d_P, mu_THETA

  term1 <- (1-rho)*P/(2*d_P)

  term2 <- a*(I_H/H)*b_V*S_V

  term3 <- mu_THETA*S_V

  return(term1 - term2 - term3)
}

calc_E_V <- function(){
  # inputs
  #   a, I_H, H, b_V, S_V, tau, mu_THETA, E_V
  #   time delayed version of I_V and S_V: tdI_H, tdS_V

  term1 <- a*(I_H/H)*b_V*S_V

  term3 <- mu_THETA*E_V

  term2 <- a*(tdI_H/H)*b_V*tdS_V*exp(-mu_THETA*tau)

  return(term1-term2-term3)
}

calc_I_V <- function(){
  # inputs
  #   a, H, b_V, mu_THETA, tau, I_V
  #   time delayed version of I_H and S_V: tdI_H, tdS_V

  term1 <- a*(tdI_H/H)*b_V*tdS_V*exp(-mu_THETA*tau)

  term2 <- mu_THETA*I_V

  return(term1-term2)
}

###############################################################################
# Run Life Cycle Model
###############################################################################
#############################
# Setup inputs for DDE
#############################
# Things that should change over time (ie, input is a vector)
#   chi_ITN
#   chi_IRS
#   I_H

time <-
chi_ITN <-
chi_IRS <-
I_H <-

## setup initial conditions
omega <- calc_omega(Parameters1 = Parameters1)
state <- calc_MosquitoEquilibrium(Parameters1 = Parameters1)
K <- calc_K(Parameters1 = Parameters1, MosEQ = state)

p <- c(Parameters1, Parameters2)
p["chi_ITN"] <- chi_ITN
p["chi_IRS"] <- chi_IRS
p["I_H"] <- I_H

#############################
# DDe Model
#############################
MosquitoLives <- function(t,state,p){
  with(as.list(c(state,p)),{

    # calculate Pr, Pf, mu_THETA
    Pr <- calc_Pr()
    Pf <- calc_Pf()
    mu_THETA <- calc_muFTheta()

    # retrieve time-delayed parts
    if(t >= tau){
      tdI_H <- I_H[t-tau+1]
      tdS_V <- lagvalue(t-tau, 4) # check this!!! idk what it might return
    } else {
      tdI_H <- 0
      tdS_V <- 0
    }

    # aquatic equations
    dE <- calc_dE()
    dL <- calc_dL()
    dP <- calc_dP()

    # adult equations
    dS_V <- calc_dS_V()
    dE_V <- calc_E_V()
    dI_V <- calc_I_V()

    # return
    return(c(dE,dL,dP,dS_V,dE_V,dI_V))
  })
}

#############################
# Solve and Plot Model
#############################
## solve mosquito life model
ML_out <- deSolve::dede(y = state, times = time, func = MosquitoLives,
                        parms = p, method = "lsoda",
                        control = list(mxhist=length(time)))

## plot using this beautiful plotting example
#?
#?
#?
#?




















