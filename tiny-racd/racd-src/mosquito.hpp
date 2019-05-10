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
#include <Rmath.h>

#include <vector>
#include <memory>

// // output: states
// extern std::vector<size_t> state_S;
// extern std::vector<size_t> state_E;
// extern std::vector<size_t> state_T;
// extern std::vector<size_t> state_D;
// extern std::vector<size_t> state_A;
// extern std::vector<size_t> state_U;
// extern std::vector<size_t> state_P;
//
// // output: pop sizes
// extern std::vector<size_t> num_All;
// extern std::vector<size_t> num_2_10;
// extern std::vector<size_t> num_0_5;
// extern std::vector<size_t> num_5_10;
// extern std::vector<size_t> num_10_15;
// extern std::vector<size_t> num_15Plus;
//
// // output: clinical incidence
// extern std::vector<size_t> cinc_All;
// extern std::vector<size_t> cinc_2_10;
// extern std::vector<size_t> cinc_0_5;
// extern std::vector<size_t> cinc_5_10;
// extern std::vector<size_t> cinc_10_15;
// extern std::vector<size_t> cinc_15Plus;
//
// // output for mosquitos
// extern std::vector<size_t> mosy_S;
// extern std::vector<size_t> mosy_E;
// extern std::vector<size_t> mosy_I;
//
// // parameters
// extern std::unordered_map<std::string,double> parameters;
//
// // current simulation time
// extern size_t tnow;
//
// // id for people
// extern size_t global_hid;
//
// // transmission: mosy -> human
//
// // psi (biting weights) and EIR for each house
// extern std::vector<double>   psi;
// extern std::vector<double>   EIR;
//
// // transmission: human -> mosy
//
// extern double      CC; // P(bite will cause infection in mosquito) --- expectation of this prob over all landscape/individual heterogeneities
// extern double      WW; // avg probability to bite and survive
// extern double      ZZ; // avg probability to bite
//
// extern void reset_states(const size_t tmax);
// extern void reset_pop(const size_t tmax);
// extern void reset_cinc(const size_t tmax);
// extern void reset_mosy(const size_t tmax);
// extern void reset_globals(const size_t tmax);

// the mosquito population struct
typedef struct mosquitos {

  /* new eggs are generated from a conditionally independent Poisson process */
  int                      EL_new;

  /* probabilities & transitions for early-stage instars ("EL","D","LL") */
  std::vector<double>         EL_probs;
  std::vector<int>         EL_transitions;

  /* probabilities & transitions for late-stage instars ("LL","D","PL") */
  std::vector<double>         LL_probs;
  std::vector<int>         LL_transitions;

  /* probabilities & transitions for pupae ("PL","D","SV_F","SV_M") */
  std::vector<double>         PL_probs;
  std::vector<int>         PL_transitions;

  /* probabilities & transitions for susceptible vectors ("SV","D","EV") */
  std::vector<double>         SV_probs;
  std::vector<int>         SV_transitions;

  /* probabilities & transitions for incubating vectors ("EV","D","IV") */
  std::vector<double>         EV_probs;
  std::vector<int>         EV_transitions;

  /* probabilities & transitions for infectious vectors ("IV","D") */
  std::vector<double>         IV_probs;
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
  mosquitos(const int EL_, const int LL_, const int PL_, const int SV_, const int EV_, const int IV_, const double K_);
  ~mosquitos();

} mosquitos;

using mosquito_ptr = std::unique_ptr<mosquitos>;


/* ################################################################################
#   daily updates
################################################################################ */

void feeding_cycle(mosquito_ptr& mosy);

void euler_step(mosquito_ptr& mosy);

#endif
