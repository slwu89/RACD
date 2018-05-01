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

#ifndef RACD_HUMAN_
#define RACD_HUMAN_

#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <functional>
#include <math.h> /* for update_lambda */

using namespace std::placeholders;

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
  void                            one_day(const int tNow);

private:

  /* hash table of pointers to member functions */
  std::unordered_map<std::string,std::function<void(const int)> >   compartment_funs;

  /* private methods */
  void                            mortality(const int tNow);      /* mortality */
  void                            S_compartment(const int tNow);  /* S: susceptible */
  void                            E_compartment(const int tNow);  /* E: latent liver-stage infection */
  void                            T_compartment(const int tNow);  /* T: treated clinical disease */
  void                            D_compartment(const int tNow);  /* D: untreated clinical disease */
  void                            A_compartment(const int tNow);  /* A: asymptomatic patent (detectable by microscopy) infection */
  void                            U_compartment(const int tNow);  /* U: asymptomatic sub-patent (not detectable by microscopy) infection */
  void                            P_compartment(const int tNow);  /* P: protection due to chemoprophylaxis treatment */
  void                            ageing();                       /* ageing */
  void                            update_immunity();              /* immunity */
  void                            update_lambda();                /* lambda: force of infection */
  void                            update_phi();                   /* phi */
  void                            update_q();                     /* q (microscopy) */

  /* parameters */
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
