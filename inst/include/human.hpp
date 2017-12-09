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

#include <Rcpp.h>
#include <string>

// typedefs and forward declarations
class house;



class human {

public:

  human();
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

  double                          epsilon; /* force of infection */
  double                          lambda; /* EIR */
  double                          phi;

  house*                          myHouse; /* raw pointer ok because house lifespan > human lifespan in memory */
};
