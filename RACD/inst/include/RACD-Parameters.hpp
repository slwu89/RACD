/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  December 2017
 #
 #  Parameters Singleton Class Definition
*/

/* ######################################################################
 # includes and foward declarations
###################################################################### */

#ifndef RACD_PARAMETERS
#define RACD_PARAMETERS

/* C++ includes */
#include <iostream>
#include <memory>

// #include "DEBUG.hpp"

/* ######################################################################
 # class declaration
###################################################################### */

class parameters {

public:

  /* constructor & destructor */
  parameters(
      double epsilon0_,
      double fT_,
      int dE_,
      int dT_,
      int dD_,
      int dA_,
      int dU_,
      int dP_,
      double cD_,
      double cT_,
      double cU_,
      double gammaI_,
      double rho_,
      double a0_,
      double sigma2_,
      double d1_,
      double dID_,
      double ID0_,
      double kappaD_,
      double uD_,
      double aD_,
      double fD0_,
      double gammaD_,
      double alphaA_,
      double alphaU_,
      double b0_,
      double b1_,
      double dB_,
      double IB0_,
      double kappaB_,
      double uB_,
      double phi0_,
      double phi1_,
      double dC_,
      double IC0_,
      double kappaC_,
      double uC_,
      double PM_,
      double dM_,
      double rW_,
      double rP_,
      double meanAge_,
      int N_);
  ~parameters();

  /* delete all copy & move semantics */
  parameters(const parameters&) = delete;
  parameters& operator=(const parameters&) = delete;
  parameters(parameters&&) = delete;
  parameters& operator=(parameters&&) = delete;

  /* accessors */
  double get_epsilon0(){ return epsilon0; };
  double get_fT(){ return fT; };
  int get_dE(){ return dE; };
  int get_dT(){ return dT; };
  int get_dD(){ return dD; };
  int get_dA(){ return dA; };
  int get_dU(){ return dU; };
  int get_dP(){ return dP; };
  double get_cD(){ return cD; };
  double get_cT(){ return cT; };
  double get_cU(){ return cU; };
  double get_gammaI(){ return gammaI; };
  double get_rho(){ return rho; };
  double get_a0(){ return a0; };
  double get_sigma2(){ return sigma2; };
  double get_d1(){ return d1; };
  double get_dID(){ return dID; };
  double get_ID0(){ return ID0; };
  double get_kappaD(){ return kappaD; };
  double get_uD(){ return uD; };
  double get_aD(){ return aD; };
  double get_fD0(){ return fD0; };
  double get_gammaD(){ return gammaD; };
  double get_alphaA(){ return alphaA; };
  double get_alphaU(){ return alphaU; };
  double get_b0(){ return b0;};
  double get_b1(){ return b1;};
  double get_dB(){ return dB;};
  double get_IB0(){ return IB0;};
  double get_kappaB(){ return kappaB;};
  double get_uB(){ return uB;};
  double get_phi0(){ return phi0;};
  double get_phi1(){ return phi1;};
  double get_dC(){ return dC;};
  double get_IC0(){ return IC0;};
  double get_kappaC(){ return kappaC;};
  double get_uC(){ return uC;};
  double get_PM(){ return PM;};
  double get_dM(){ return dM;};
  double get_rW(){ return rW;};
  double get_rP(){ return rP;};
  double get_meanAge(){ return meanAge;};
  int get_N(){ return N;};
  double get_mu(){ return mu; };

private:
  /* parameters */
  double                      epsilon0; /* Mean EIR for adults (per day) */
  double                      fT; /* Proportion of clinical disease cases successfully treated */

  /* Model parameters taken from Griffin et al. (2014) */
  /* Human infection durations */
  int                         dE; /* Duration of latent period (days) */
  int                         dT; /* Duration of treated clinical disease (days) */
  int                         dD; /* Duration of untreated clinical disease (days) */
  int                         dA; /* Duration of patent infection (days) */
  int                         dU; /* Duration of sub-patent infection (days) (fitted) */
  int                         dP; /* Duration of prophylactic protection following treatment (days) */

  /* Infectiousness of humans to mosquitoes */
  double                      cD; /* Infectiousness with untreated disease & no immunity (fitted) */
  double                      cT; /* Infectiousness after treatment */
  double                      cU; /* Infectiousness with sub-patent infection (fitted) */
  double                      gammaI; /* Relates infectiousness to probability of detection (fitted) */

  /* Age and heterogeneity parameters */
  double                      rho; /* Age-dependent biting parameter */
  double                      a0; /* Age-dependent biting parameter (years) */
  double                      sigma2; /* Variance of log of heterogeneity in biting rates */

  /* Effect of immunity on reducing probability of detection */
  double                      d1; /* Probability of detection with maximum immunity (fitted) */
  double                      dID; /* Inverse of decay rate (days) */
  double                      ID0; /* Immunity scale parameter (fitted) */
  double                      kappaD; /* Immunity shape parameter (fitted) */
  double                      uD; /* Duration in which immunity is not boosted (days) (fitted) */
  double                      aD; /* Scale parameter relating age to immunity (years) (fitted) */
  double                      fD0; /* Parameter relating age to immunity (fitted) */
  double                      gammaD; /* Shape parameter relating age to immunity (fitted) */
  double                      alphaA; /* PCR prevalence parameter (fitted) */
  double                      alphaU; /* PCR prevalence parameter (fitted) */

  /* Immunity reducing probability of infection */
  double                      b0; /* Probabiliy with no immunity (fitted) */
  double                      b1; /* Maximum relative reduction */
  double                      dB; /* Inverse of decay rate (days) */
  double                      IB0; /* Scale parameter (fitted) */
  double                      kappaB; /* Shape parameter (fitted) */
  double                      uB; /* Duration in which immunity is not boosted (days) (fitted) */

  /* Immunity reducing probability of clinical disease */
  double                      phi0; /* Probability with no immunity */
  double                      phi1; /* Maximum relative reduction */
  double                      dC; /* Inverse decay rate (days) */
  double                      IC0; /* Scale parameter */
  double                      kappaC; /* Shape parameter */
  double                      uC; /* Duration in which immunity is not boosted (days) */
  double                      PM; /* New-born immunity relative to mother's immunity */
  double                      dM; /* Inverse decay rate of maternal immunity (days) */

  /* Case detection (recorded incidence relative to daily active case detection) */
  double                      rW; /* Weekly active case detection */
  double                      rP; /* Weekly passive case detection */

  /* Demographic parameters */
  double                      meanAge; /* Mean age in Tanzania (males and females, years) */
  int                         N; /* Village population size */
  double                      mu; /* Daily death rate as a function of mean age in years */
};

#endif
