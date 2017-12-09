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

#include <Rcpp.h>
#include <iostream>

#include "debug.hpp"

class RACD_Parameters {
public:


private:

  /* constructor */

  /* destructor */


  /* parameters */
  double              epsilon0; /* Mean EIR for adults (per day) */
  double              fT; /* Proportion of clinical disease cases successfully treated */

  /* Model parameters taken from Griffin et al. (2014) */
  /* Human infection durations */
  int                 dE; /* Duration of latent period (days) */
  int                 dT; /* Duration of treated clinical disease (days) */
  int                 dD; /* Duration of untreated clinical disease (days) */
  int                 dA; /* Duration of patent infection (days) */
  int                 dU; /* Duration of sub-patent infection (days) (fitted) */
  int                 dP; /* Duration of prophylactic protection following treatment (days) */

  /* Infectiousness of humans to mosquitoes */
  double              cD; /* Infectiousness with untreated disease & no immunity (fitted) */
  double              cT; /* Infectiousness after treatment */
  double              cU; /* Infectiousness with sub-patent infection (fitted) */
  double              gammaI; /* Relates infectiousness to probability of detection (fitted) */

  /* Age and heterogeneity parameters */
  double              rho; /* Age-dependent biting parameter */
  double              a0; /* Age-dependent biting parameter (years) */
  double              sigma2; /* Variance of log of heterogeneity in biting rates */

  /* Effect of immunity on reducing probability of detection */
  double              d1; /* Probability of detection with maximum immunity (fitted) */
  double              dID; /* Inverse of decay rate (days) */
  double              ID0; /* Immunity scale parameter (fitted) */
  double              kappaD; /* Immunity shape parameter (fitted) */
  double              uD; /* Duration in which immunity is not boosted (days) (fitted) */
  double              aD; /* Scale parameter relating age to immunity (years) (fitted) */
  double              fD0; /* Parameter relating age to immunity (fitted) */
  double              gammaD; /* Shape parameter relating age to immunity (fitted) */
  double              alphaA; /* PCR prevalence parameter (fitted) */
  double              alphaU; /* PCR prevalence parameter (fitted) */

  /* Immunity reducing probability of infection */
  double              b0; /* Probabiliy with no immunity (fitted) */
  double              b1; /* Maximum relative reduction */
  double              dB; /* Inverse of decay rate (days) */
  double              IB0; /* Scale parameter (fitted) */
  double              kappaB; /* Shape parameter (fitted) */
  double              uB; /* Duration in which immunity is not boosted (days) (fitted) */

  /* Immunity reducing probability of clinical disease */
  double              phi0; /* Probability with no immunity */
  double              phi1; /* Maximum relative reduction */
  double              dC /* Inverse decay rate (days) */
  double              IC0; /* Scale parameter */
  double              kappaC; /* Shape parameter */
  double              uC; /* Duration in which immunity is not boosted (days) */
  double              PM; /* New-born immunity relative to mother's immunity */
  double              dM; /* Inverse decay rate of maternal immunity (days) */

  /* Case detection (recorded incidence relative to daily active case detection) */
  double              rW; /* Weekly active case detection */
  double              rP; /* Weekly passive case detection */

  /* Demographic parameters */
  double              meanAge; /* Mean age in Tanzania (males and females, years) */
  int                 N; /* Village population size */

  /* Geographic parameters */
  double              meanNumPeoplePerHouse; /* Mean number of people per house (from Misungu data set) */
  int                 numHousesPerBreedingSite; /* Number of houses per breeding site */


};





class GlobalClass
{
    int m_value;
    static GlobalClass *s_instance;
    GlobalClass(int v = 0){
        std::cout << "global singleton being born at " << this << std::endl;
        m_value = v;
    };
    ~GlobalClass(){
        std::cout << "global singleton being killed at " << this << std::endl;
    };
public:
    int get_value()
    {
        return m_value;
    }
    void set_value(int v)
    {
        m_value = v;
    }
    void kill_me(){
        std::cout << "global singleton is dying!" << std::endl;
        delete this;
    }
    static GlobalClass *instance()
    {
        if (!s_instance)
            s_instance = new GlobalClass;
        return s_instance;
    }
};



friend std::ostream& operator<<(std::ostream& os, tentPAR& tent) {
  return os << "tentPAR: " << " t0: " << tent.t0 << " pfid: " << tent.pfid << " gr: " << tent.gr << " dr: " << tent.dr << "\n" <<
    " MZ0 " << tent.MZ0 << " peakD " << tent.peakD << " mxPD: " << tent.mxPD << " tEnd: " << tent.tEnd << std::endl;
}
