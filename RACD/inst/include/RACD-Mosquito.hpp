/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Marshall Lab
 #  February 2019
 #
 #  Mosquito Habitat Class Definition
*/

#ifndef RACD_HABITAT
#define RACD_HABITAT

#include <Rcpp.h>

#include <iostream>
#include <vector>

#include <math.h>
#include <Rmath.h> // for rmultinom

/* alias and forward declarations */
class village;

/* mosquito */
class mosquito_habitat {
public:

  /* constructor & destructor */
  mosquito_habitat(const int EL_, const int LL_, const int PL_, const int SV_, const int EV_, const int IV_, const double K_,
    const std::vector<double>& psiR, village* const village_ptr_);
  ~mosquito_habitat();

  /* delete all copy semantics */
  mosquito_habitat(const mosquito_habitat&) = delete;
  mosquito_habitat& operator=(const mosquito_habitat&) = delete;

  /* default move semantics */
  mosquito_habitat(mosquito_habitat&&);
  mosquito_habitat& operator=(mosquito_habitat&&);

  /* simulation */
  void feeding_cycle();
  void euler_step(const double tnow);

  /* accessors */
  int get_EL(){return EL;}
  int get_LL(){return LL;}
  int get_PL(){return PL;}
  int get_SV(){return SV;}
  int get_EV(){return EV;}
  int get_IV(){return IV;}

private:

  /* biting */
  std::vector<double>         psiWeight; /* weights (assuming all houses have people, should still sum to 1) */
  std::vector<double>         psi; /* vector of biting distribution on houses */
  void                        normalize_psi();

  village* const              village_ptr; /* raw pointer ok because house lifespan > human lifespan in memory */

  /* new eggs are generated from a conditionally independent Poisson process */
  int                         EL_new;

  /* probabilities & transitions for early-stage instars ("EL","D","LL") */
  std::vector<double>         EL_probs;
  std::vector<int>            EL_transitions;

  /* probabilities & transitions for late-stage instars ("LL","D","PL") */
  std::vector<double>         LL_probs;
  std::vector<int>            LL_transitions;

  /* probabilities & transitions for pupae ("PL","D","SV_F","SV_M") */
  std::vector<double>         PL_probs;
  std::vector<int>            PL_transitions;

  /* probabilities & transitions for susceptible vectors ("SV","D","EV") */
  std::vector<double>         SV_probs;
  std::vector<int>            SV_transitions;

  /* probabilities & transitions for incubating vectors ("EV","D","IV") */
  std::vector<double>         EV_probs;
  std::vector<int>            EV_transitions;

  /* probabilities & transitions for infectious vectors ("IV","D") */
  std::vector<double>         IV_probs;
  std::vector<int>            IV_transitions;

  /* state space */
  int                         EL;
  int                         LL;
  int                         PL;
  int                         SV;
  int                         EV;
  int                         IV;

  /* carrying capacity */
  double                      K;

  /* outbound bites */
  std::vector<int>            EIR_out;

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

};


#endif
