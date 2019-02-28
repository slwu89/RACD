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
#include <unordered_map>

#include <math.h>

/* alias and forward declarations */
class village;


/* ################################################################################
 * template class: int for stochastic, double for deterministic
################################################################################ */

template <typename T>
class mosquito_habitat {
public:

  /* constructor & destructor */
  mosquito_habitat(const T EL_, const T LL_, const T PL_, const T SV_, const T EV_, const T IV_, const double K_, const Rcpp::List pars_);
  ~mosquito_habitat();

  /* delete all copy semantics */
  mosquito_habitat(const mosquito_habitat&) = delete;
  mosquito_habitat& operator=(const mosquito_habitat&) = delete;

  /* default move semantics */
  mosquito_habitat(mosquito_habitat&&);
  mosquito_habitat& operator=(mosquito_habitat&&);

  /* simulation */
  void euler_step(const double tnow, const double dt);

  /* accessors */
  T get_EL(){return EL;}
  T get_LL(){return LL;}
  T get_PL(){return PL;}
  T get_SV(){return SV;}
  T get_EV(){return EV;}
  T get_IV(){return IV;}

private:

  /* default data members */
  size_t                      habitatID;
  std::vector<double>         psi; /* vector of biting distribution on houses */
  village*                    village_ptr; /* raw pointer ok because house lifespan > human lifespan in memory */

  /* new eggs are generated from a conditionally independent Poisson process */
  T                     EL_new;

  /* probabilities & transitions for early-stage instars ("EL","D","LL") */
  std::vector<double>   EL_probs;
  std::vector<T>        EL_transitions;

  /* probabilities & transitions for late-stage instars ("LL","D","PL") */
  std::vector<double>   LL_probs;
  std::vector<T>        LL_transitions;

  /* probabilities & transitions for pupae ("PL","D","SV_F","SV_M") */
  std::vector<double>   PL_probs;
  std::vector<T>        PL_transitions;

  /* probabilities & transitions for susceptible vectors ("SV","D","EV") */
  std::vector<double>   SV_probs;
  std::vector<T>        SV_transitions;

  /* probabilities & transitions for incubating vectors ("EV","D","IV") */
  std::vector<double>   EV_probs;
  std::vector<T>        EV_transitions;

  /* probabilities & transitions for infectious vectors ("IV","D") */
  std::vector<double>   IV_probs;
  std::vector<T>        IV_transitions;

  /* state space */
  T   EL;
  T   LL;
  T   PL;
  T   SV;
  T   EV;
  T   IV;

  /* carrying capacity */
  double K;

  /* parameters */
  std::unordered_map<std::string, double> pars;

};


#endif
