/* ################################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John M. Marshall
#   December 2017
#
#   Human Class Definition
################################################################################ */

#ifndef RACD_HUMAN_
#define RACD_HUMAN_

#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <functional>
#include <math.h> /* for update_lambda */

using namespace std::placeholders;

// #include "DEBUG.hpp"

/* alias and forward declarations */
class house;

/* human */
class human {

public:

  /* constructor */
  human(const int humanID_,
        const double age_,
        const bool alive_,
        const std::string& state_,
        const int daysLatent_,
        const double IB_,
        const double ID_,
        const double ICA_,
        const double ICM_,
        const double bitingHet_,
        const double epsilon_,
        const double lambda_,
        const double phi_,
        const double prDetectAMic_,
        const double prDetectAPCR_,
        const double prDetectUPCR_,
        house* house_ptr_
      );

  /* destructor */
  ~human();

  /* accessors */
  double                          get_age(){ return age; };
  double                          get_ICA(){ return ICA; };
  bool const&                     get_alive() const { return alive; };
  double                          get_c(){ return c; };

  /* Simulation Methods */

  /* daily simulation */
  void                            one_day(const int tNow);

private:

  /* hash table of pointers to member functions */
  std::unordered_map<std::string,std::function<void(const int)> >   compartment_funs;
  std::unordered_map<std::string,std::function<void()> >            infectiousness_funs;

  /* private methods */
  void                            S_compartment(const int tNow);  /* S: susceptible */
  void                            E_compartment(const int tNow);  /* E: latent liver-stage infection */
  void                            T_compartment(const int tNow);  /* T: treated clinical disease */
  void                            D_compartment(const int tNow);  /* D: untreated clinical disease */
  void                            A_compartment(const int tNow);  /* A: asymptomatic patent (detectable by microscopy) infection */
  void                            U_compartment(const int tNow);  /* U: asymptomatic sub-patent (not detectable by microscopy) infection */
  void                            P_compartment(const int tNow);  /* P: protection due to chemoprophylaxis treatment */

  void                            mortality(const int tNow);      /* mortality */
  void                            ageing();                       /* ageing */
  void                            update_immunity();              /* immunity */
  void                            update_lambda();                /* lambda: force of infection */
  void                            update_phi();                   /* phi */
  void                            update_q();                     /* q (microscopy) */

  void                            infectiousness_S();             /* infectiousness to mosquitoes */
  void                            infectiousness_E();
  void                            infectiousness_T();
  void                            infectiousness_D();
  void                            infectiousness_A();
  void                            infectiousness_U();
  void                            infectiousness_P();

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
  double                          c; /* infectiousness to mosquitoes */

  double                          prDetectAMic; /* q: probability of asymptomatic infection detected by microscopy */
  double                          prDetectAPCR; /* probability of detection by PCR for A (patent) asymptomatic */
  double                          prDetectUPCR; /* probability of detection by PCR for U (subpatent) asymptomatic */

  house*                          house_ptr; /* raw pointer ok because house lifespan > human lifespan in memory */
};

#endif
