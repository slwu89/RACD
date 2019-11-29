/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  November 2019
 #
 #  Mosquitoes
*/

#ifndef MOSQUITO_HPP
#define MOSQUITO_HPP

#include <vector>
#include <memory>
#include <string>
#include <unordered_map>

#include <Rcpp.h>
#include <Rmath.h>


/* ################################################################################
#   the mosquito population struct
################################################################################ */

typedef struct mosquitos {

  /* new eggs are generated from a conditionally independent Poisson process */
  int                      EL_new;

  /* probabilities & transitions for early-stage instars ("EL","D","LL") */
  std::vector<double>      EL_probs;
  std::vector<int>         EL_transitions;

  /* probabilities & transitions for late-stage instars ("LL","D","PL") */
  std::vector<double>      LL_probs;
  std::vector<int>         LL_transitions;

  /* probabilities & transitions for pupae ("PL","D","SV_F","SV_M") */
  std::vector<double>      PL_probs;
  std::vector<int>         PL_transitions;

  /* probabilities & transitions for susceptible vectors ("SV","D","EV") */
  std::vector<double>      SV_probs;
  std::vector<int>         SV_transitions;

  /* probabilities & transitions for incubating vectors ("EV","D","IV") */
  std::vector<double>      EV_probs;
  std::vector<int>         EV_transitions;

  /* probabilities & transitions for infectious vectors ("IV","D") */
  std::vector<double>      IV_probs;
  std::vector<int>         IV_transitions;

  /* state space */
  int                      EL;
  int                      LL;
  int                      PL;
  int                      SV;
  int                      EV;
  int                      IV;

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
  mosquitos(const int EL_, const int LL_, const int PL_, const int SV_, const int EV_, const int IV_, const double K_,
            const double lambdaV_);
  ~mosquitos();

} mosquitos;


/* ################################################################################
#   initialize mosquito and return to R
################################################################################ */

Rcpp::XPtr<mosquitos> init_mosquitos(
  const int EL_,
  const int LL_,
  const int PL_,
  const int SV_,
  const int EV_,
  const int IV_,
  const double K_,
  const double lambdaV_
);


/* ################################################################################
#   daily updates
################################################################################ */

// feeding cycle needs to return the vector of EIR for houses;
// bite_probs: vector of WW,ZZ,CC
// psi: vector of probabilities to distribute bites to houses
std::vector<double> feeding_cycle(SEXP mosy,
  const std::vector<double>& bite_probs,
  const std::vector<double>& psi,
  const Rcpp::NumericVector& parameters
);

void euler_step(SEXP mosy,
  const Rcpp::NumericVector& parameters
);


/* ################################################################################
#   track output
################################################################################ */

std::vector<int> track_mosquito(SEXP mosy);

double           track_lambdaV(SEXP mosy);

Rcpp::List       mosquito_2list(SEXP mosy);

#endif
