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

/* typedefs and forward declarations */
class house;

/* human */
class human {

public:

  /* constructor */
  human();

  /* destructor */
  ~human();

private:

  double                          age;
  bool                            alive;
  std::string                     state; /* S,T,D,A,U,P */
  int                             daysLatent; /* days of latent infection */


  double                          IB; /* pre-erythrocytic immunity (reduces probability of infection following infectious challenge) */
  double                          ID; /* detection immunity (blood-stage immunity, reduces the probability of detection and reduces infectiousness to mosquitoes) */
  double                          ICA; /* acquired clinical immunity (reduces the probability of clinical disease, acquired from previous exposure) */
  double                          ICM; /* maternal clinical immunity (reduces the probability of clinical disease, acquired maternally) */

  double                          bitingHet; /* zeta: baseline biting heterogeneity */
  double                          epsilon; /* force of infection */
  double                          lambda; /* EIR */
  double                          phi; /* probability of acquiring clinical disease */

  double                          prDetectAMic; /* q: probability of asymptomatic infection detected by microscopy */
  double                          prDetectAPCR; /* probability of detection by PCR for A (patent) asymptomatic */
  double                          prDetectUPCR; /* probability of detection by PCR for U (subpatent) asymptomatic */

  house*                          myHouse; /* raw pointer ok because house lifespan > human lifespan in memory */
};

#endif
