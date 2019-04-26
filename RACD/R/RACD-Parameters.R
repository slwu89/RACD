###############################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John Marshall
#   December 2017
#   RACD Model Parameters
#
###############################################################################

#' Generate Parameters for RACD Model
#'
#' Return a named vector of parameters.
#'
#' @param epsilon0 Mean EIR for adults (per day)
#' @param fT Proportion of clinical disease cases successfully treated
#' @param dE Duration of latent period (days)
#' @param dT Duration of treated clinical disease (days)
#' @param dD Duration of untreated clinical disease (days)
#' @param dA Duration of patent infection (days)
#' @param dU Duration of sub-patent infection (days) (fitted)
#' @param dP Duration of prophylactic protection following treatment (days)
#' @param cD Infectiousness with untreated disease & no immunity (fitted)
#' @param cT Infectiousness after treatment
#' @param cU Infectiousness with sub-patent infection (fitted)
#' @param gammaI Relates infectiousness to probability of detection (fitted)
#' @param rho Age-dependent biting parameter
#' @param a0 Age-dependent biting parameter (years)
#' @param sigma2 Variance of log of heterogeneity in biting rates
#' @param d1 Probability of detection with maximum immunity (fitted)
#' @param dID Inverse of decay rate (days)
#' @param ID0 Immunity scale parameter (fitted)
#' @param kappaD Immunity shape parameter (fitted)
#' @param uD Duration in which immunity is not boosted (days) (fitted)
#' @param aD Scale parameter relating age to immunity (years) (fitted)
#' @param fD0 Parameter relating age to immunity (fitted)
#' @param gammaD Shape parameter relating age to immunity (fitted)
#' @param alphaA PCR prevalence parameter (fitted)
#' @param alphaU PCR prevalence parameter (fitted)
#' @param b0 Probabiliy with no immunity (fitted)
#' @param b1 Maximum relative reduction
#' @param dB Inverse of decay rate (days)
#' @param IB0 Scale parameter (fitted)
#' @param kappaB Shape parameter (fitted)
#' @param uB Duration in which immunity is not boosted (days) (fitted)
#' @param phi0 Probability with no immunity
#' @param phi1 Maximum relative reduction
#' @param dC Inverse decay rate (days)
#' @param IC0 Scale parameter
#' @param kappaC Shape parameter
#' @param uC Duration in which immunity is not boosted (days)
#' @param PM New-born immunity relative to mother's immunity
#' @param dM Inverse decay rate of maternal immunity (days)
#' @param rW Weekly active case detection
#' @param rP Weekly passive case detection
#' @param meanAge Mean age in Tanzania (males and females, years)
#' @param N Village population size
#' @param meanNumPeoplePerHouse Mean number of people per house (from Misungu data set)
#' @param numHousesPerBreedingSite Number of houses per breeding site
#'
#' @export
RACD_Parameters <- function(
  ## Variable model parameters:
  epsilon0 = 5/365, # Mean EIR for adults (per day)
  fT = 0.4, # Proportion of clinical disease cases successfully treated

  ## Model parameters taken from Griffin et al. (2014):
  ## Human infection durations:
  dE = 12, # Duration of latent period (days)
  dT = 5, # Duration of treated clinical disease (days)
  dD = 5, # Duration of untreated clinical disease (days)
  dA = 200, # Duration of patent infection (days)
  dU = 110, # Duration of sub-patent infection (days) (fitted)
  dP = 25, # Duration of prophylactic protection following treatment (days)

  ## Infectiousness of humans to mosquitoes:
  cD = 0.068, # Infectiousness with untreated disease & no immunity (fitted)
  cT = 0.32 * 0.068, # Infectiousness after treatment
  cU = 0.0062, # Infectiousness with sub-patent infection (fitted)
  gammaI = 1.82, # Relates infectiousness to probability of detection (fitted)

  ## Age and heterogeneity parameters:
  rho = 0.85, # Age-dependent biting parameter
  a0 = 8, # Age-dependent biting parameter (years)
  sigma2 = 1.67, # Variance of log of heterogeneity in biting rates

  ## Effect of immunity on reducing probability of detection:
  d1 = 0.161, # Probability of detection with maximum immunity (fitted)
  dID = 10*365, # Inverse of decay rate (days)
  ID0 = 1.58, # Immunity scale parameter (fitted)
  kappaD = 0.477, # Immunity shape parameter (fitted)
  uD = 9.45, # Duration in which immunity is not boosted (days) (fitted)
  aD = 21.9, # Scale parameter relating age to immunity (years) (fitted)
  fD0 = 0.0071, # Parameter relating age to immunity (fitted)
  gammaD = 4.81, # Shape parameter relating age to immunity (fitted)
  alphaA = 0.757, # PCR prevalence parameter (fitted)
  alphaU = 0.186, # PCR prevalence parameter (fitted)

  ## Immunity reducing probability of infection:
  b0 = 0.590, # Probabiliy with no immunity (fitted)
  b1 = 0.5, # Maximum relative reduction
  dB = 10*365, # Inverse of decay rate (days)
  IB0 = 43.9, # Scale parameter (fitted)
  kappaB = 2.16, # Shape parameter (fitted)
  uB = 7.20, # Duration in which immunity is not boosted (days) (fitted)

  ## Immunity reducing probability of clinical disease:
  phi0 = 0.792, # Probability with no immunity
  phi1 = 0.00074, # Maximum relative reduction
  dC = 30*365, # Inverse decay rate (days)
  IC0 = 18.0, # Scale parameter
  kappaC = 2.37, # Shape parameter
  uC = 6.06, # Duration in which immunity is not boosted (days)
  PM = 0.774, # New-born immunity relative to mother's immunity
  dM = 67.7, # Inverse decay rate of maternal immunity (days)

  ## Case detection (recorded incidence relative to daily active case
  ## detection):
  rW = 0.723, # Weekly active case detection
  rP = 0.342, # Weekly passive case detection

  ## Demographic parameters:
  meanAge = 17.4, # Mean age in Tanzania (males and females, years)
  N = 200, # Village population size

  ## Geographic parameters:
  meanNumPeoplePerHouse = 6.5, # Mean number of people per house (from Misungu data set)
  numHousesPerBreedingSite = 5,

  # ENTO PARS
  ## Mosquito life cycle parameters:
  beta = 21.19, # Number of eggs laid per day by female mosquito
  muEL = 0.034, # Early larval instar daily mortality
  muLL = 0.035, # Late larval instar daily mortality
  muPL = 0.25, # Pupal daily mortality
  durEL = 6.64, # Duration of early instar stage
  durLL = 3.72, # Duration of late instar stage
  durPL = 0.64, # Duration of pupal stage
  durEV = 10, # Duration of latent period in mosquito (days)
  gamma = 13.25, # Effect of density-dependence on late instars relative to early instars
  tau1 = 0.68, # Time spent foraging for a blood meal at 0% ITN coverage
  tau2 = 2.32, # Time spent resting and ovipositing by a mosquito

  muV = 1/7.6, # Adult mosquito daily mortality
  Q0 = 0.92, # Human blood index
  phiB = 0.89, # Proportion of bites on a person while they are in bed
  phiI = 0.97, # Proportion of bites on a person while they are indoors
  rITN = 0.56, # Probability of mosquito repeating a feeding attempt due to IRS
  sITN = 0.03, # Probability of mosquito feeding and surviving in presence of ITNs
  rIRS = 0.60, # Probability of mosquito repeating a feeding attempt due to IRS
  sIRS = 0 # Probability of mosquito feeding and surviving in presence of IRS

){

  delta <- 1 / (tau1+tau2)
  eggOV <- beta*(exp(muV/delta)-1.0)/muV

  return(
    c(
      epsilon0 = epsilon0,
      fT = fT,
      dE = dE,
      dT = dT,
      dD = dD,
      dA = dA,
      dU = dU,
      dP = dP,
      cD = cD,
      cT = cT,
      cU = cU,
      gammaI = gammaI,
      rho = rho,
      a0 = a0,
      sigma2 = sigma2,
      d1 = d1,
      dID = dID,
      ID0 = ID0,
      kappaD = kappaD,
      uD = uD,
      aD = aD,
      fD0 = fD0,
      gammaD = gammaD,
      alphaA = alphaA,
      alphaU = alphaU,
      b0 = b0,
      b1 = b1,
      dB = dB,
      IB0 = IB0,
      kappaB = kappaB,
      uB = uB,
      phi0 = phi0,
      phi1 = phi1,
      dC = dC,
      IC0 = IC0,
      kappaC = kappaC,
      uC = uC,
      PM = PM,
      dM = dM,
      rW = rW,
      rP = rP,
      meanAge = meanAge,
      N = N,
      meanNumPeoplePerHouse = meanNumPeoplePerHouse,
      numHousesPerBreedingSite = numHousesPerBreedingSite,

      # ENTO
      beta = beta,
      muEL = muEL,
      muLL = muLL,
      muPL = muPL,
      durEL = durEL,
      durLL = durLL,
      durPL = durPL,
      durEV = durEV,
      gamma = gamma,
      tau1 = tau1,
      tau2 = tau2,

      muV = muV,
      Q0 = Q0,
      phiB = phiB,
      phiI = phiI,
      rITN = rITN,
      sITN = sITN,
      rIRS = rIRS,
      sIRS = sIRS,
      eggOV = eggOV
    )
  )
}
