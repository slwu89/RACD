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

#ifndef _RACD_PARAMETERS_
#define _RACD_PARAMETERS_

#include <Rcpp.h>
#include <iostream>
#include <memory>

#include "DEBUG.hpp"

/* typedefs and forward declarations */
class prng;
typedef std::unique_ptr<prng> prng_ptr;


class RACD_Parameters {

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

  /* Geographic parameters */
  double                      meanNumPeoplePerHouse; /* Mean number of people per house (from Misungu data set) */
  int                         numHousesPerBreedingSite; /* Number of houses per breeding site */

  /* pseudo-random number generation */
  prng_ptr                        prng_member;

  /* singleton instance */
  static RACD_Parameters      *RACD_Parameters_instance;

  /* constructor */
  RACD_Parameters();

  /* destructor */
  ~RACD_Parameters();

public:

  /* set parameter values */
  void set_values(
    double _epsilon0,
    double _fT,
    int _dE,
    int _dT,
    int _dD,
    int _dA,
    int _dU,
    int _dP,
    double _cD,
    double _cT,
    double _cU,
    double _gammaI,
    double _rho,
    double _a0,
    double _sigma2,
    double _d1,
    double _dID,
    double _ID0,
    double _kappaD,
    double _uD,
    double _aD,
    double _fD0,
    double _gammaD,
    double _alphaA,
    double _alphaU,
    double _b0,
    double _b1,
    double _dB,
    double _IB0,
    double _kappaB,
    double _uB,
    double _phi0,
    double _phi1,
    double _dC,
    double _IC0,
    double _kappaC,
    double _uC,
    double _PM,
    double _dM,
    double _rW,
    double _rP,
    double _meanAge,
    int _N,
    double _meanNumPeoplePerHouse,
    int _numHousesPerBreedingSite
  );

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
  double get_meanNumPeoplePerHouse(){ return meanNumPeoplePerHouse;};
  int get_numHousesPerBreedingSite(){ return numHousesPerBreedingSite;};

  /* pseudo-random number generation */
  void set_prng(const uint_least32_t &seed);
  prng* get_prng();

  /* suicide at end of program */
  void suicide();

  /* return instance */
  static RACD_Parameters* instance();

};

#endif
