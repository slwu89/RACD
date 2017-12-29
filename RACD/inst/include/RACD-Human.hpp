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
 #  Human Class Definition
*/

#ifndef _RACD_HUMAN_
#define _RACD_HUMAN_

#include <Rcpp.h>
#include <string>
#include <math.h> /* for update_lambda */

// #include "DEBUG.hpp"

/* typedefs and forward declarations */
class house;

/* human */
class human {

public:

  /* constructor */
  human(const int& _humanID,
        const double& _age,
        const bool& _alive,
        const std::string& _state,
        const int& _daysLatent,
        const double& _IB,
        const double& _ID,
        const double& _ICA,
        const double& _ICM,
        const double& _bitingHet,
        const double& _epsilon,
        const double& _lambda,
        const double& _phi,
        const double& _prDetectAMic,
        const double& _prDetectAPCR,
        const double& _prDetectUPCR,
        house* _house_ptr
      );

  /* destructor */
  ~human();

  /* accessors */
  double                          get_age(){ return age; };
  double                          get_ICA(){ return ICA; };
  bool const&                     get_alive() const { return alive; };

  /* Simulation Methods */

  /* daily simulation */
  void                            one_day(const int& tNow);

  /* mortality */
  void                            mortality(const int& tNow);

  /* S: susceptible */
  void                            S_compartment(const int& tNow);

  /* E: latent period */
  void                            E_compartment(const int& tNow);

  /* T: treated clinical disease */
  void                            T_compartment(const int& tNow);

  /* D: untreated clinical disease */
  void                            D_compartment(const int& tNow);

  /* A: asymptomatic patent (detectable by microscopy) infection */
  void                            A_compartment(const int& tNow);

  /* U: asymptomatic sub-patent (not detectable by microscopy) infection */
  void                            U_compartment(const int& tNow);

  /* P: protection due to chemoprophylaxis treatment */
  void                            P_compartment(const int& tNow);

  /* ageing */
  void                            ageing();

  /* immunity */
  void                            update_immunity();

  /* lambda */
  void                            update_lambda();

  /* phi */
  void                            update_phi();

  /* q (microscopy) */
  void                            update_q();

private:

  int                             humanID; /* my ID */

  double                          age; /* age (in years) */
  bool                            alive; /* alive? */
  std::string                     state; /* S,T,D,A,U,P */
  int                             daysLatent; /* days of latent infection */


  double                          IB; /* pre-erythrocytic immunity (reduces probability of infection following infectious challenge) */
  double                          ID; /* detection immunity (blood-stage immunity, reduces the probability of detection and reduces infectiousness to mosquitoes) */
  double                          ICA; /* acquired clinical immunity (reduces the probability of clinical disease, acquired from previous exposure) */
  double                          ICM; /* maternal clinical immunity (reduces the probability of clinical disease, acquired maternally) */

  double                          bitingHet; /* zeta: baseline biting heterogeneity */
  double                          epsilon; /* EIR */
  double                          lambda; /* force of infection */
  double                          phi; /* probability of acquiring clinical disease */

  double                          prDetectAMic; /* q: probability of asymptomatic infection detected by microscopy */
  double                          prDetectAPCR; /* probability of detection by PCR for A (patent) asymptomatic */
  double                          prDetectUPCR; /* probability of detection by PCR for U (subpatent) asymptomatic */

  house*                          house_ptr; /* raw pointer ok because house lifespan > human lifespan in memory */
};

#endif
