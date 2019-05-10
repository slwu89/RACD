/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  May 2019
 #
 #  the mosquitos
*/

#ifndef mosquitos_hpp
#define mosquitos_hpp

#include <Rcpp.h>

#include <vector>
#include <memory>

// the mosquito population struct
typedef struct mosquitos {

  /* new eggs are generated from a conditionally independent Poisson process */
  size_t                      EL_new;

  /* probabilities & transitions for early-stage instars ("EL","D","LL") */
  std::vector<double>         EL_probs;
  std::vector<size_t>         EL_transitions;

  /* probabilities & transitions for late-stage instars ("LL","D","PL") */
  std::vector<double>         LL_probs;
  std::vector<size_t>         LL_transitions;

  /* probabilities & transitions for pupae ("PL","D","SV_F","SV_M") */
  std::vector<double>         PL_probs;
  std::vector<size_t>         PL_transitions;

  /* probabilities & transitions for susceptible vectors ("SV","D","EV") */
  std::vector<double>         SV_probs;
  std::vector<size_t>         SV_transitions;

  /* probabilities & transitions for incubating vectors ("EV","D","IV") */
  std::vector<double>         EV_probs;
  std::vector<size_t>         EV_transitions;

  /* probabilities & transitions for infectious vectors ("IV","D") */
  std::vector<double>         IV_probs;
  std::vector<size_t>         IV_transitions;

  /* state space */
  size_t                      EL;
  size_t                      LL;
  size_t                      PL;
  size_t                      SV;
  size_t                      EV;
  size_t                      IV;

  /* carrying capacity */
  double                      K;

  /* feeding cycle vital rates */
  double                      W; // P(successful feed)
  double                      Z; // P(repelled w/out feed)
  double                      f; // feeding rate
  double                      mu; // death rate
  double                      p1; // P(survive foraging)
  double                      p2; // P(survive ovipositing)
  double                      Q; // proportion of successful bites on humans
  double                      a; // HBR
  double                      lambdaV; // FOI on mosy
  double                      beta; // eggs/day/mosquito

  /* constructor & destructor */
  mosquitos(const size_t EL_, const size_t LL_, const size_t PL_, const size_t SV_, const size_t EV_, const size_t IV_, const double K_);
  ~mosquitos();

} mosquitos;

using mosquito_ptr = std::unique_ptr<mosquitos>;


/* ################################################################################
#   daily updates
################################################################################ */

void feeding_cycle(mosquito_ptr& mosy);

void euler_step(mosquito_ptr& mosy);

#endif
