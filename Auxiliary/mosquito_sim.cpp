/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  February 2019
 #
 #  Mosquito Habitat
*/

// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>

#include <iostream>
#include <vector>
// #include <unordered_map>

#include <math.h>

/* template for int and double */
template <typename T>
class mosquito_habitat {
public:

  /* constructor & destructor */
  mosquito_habitat();
  ~mosquito_habitat();

  /* delete all copy semantics */
  mosquito_habitat(const mosquito_habitat&) = delete;
  mosquito_habitat& operator=(const mosquito_habitat&) = delete;

  /* default move semantics */
  mosquito_habitat(mosquito_habitat&&);
  mosquito_habitat& operator=(mosquito_habitat&&);

  /* simulation */
  void euler_step(const double dt);

private:

  /* new eggs are generated from a conditionally independent Poisson process */
  T   EL_new;

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

  double K;

  /* parameters */
  static std::unordered_map<std::string, double> pars;

  /* intervention parameters */

};

/* constructor */
template <typename T>
mosquito_habitat<T>::mosquito_habitat() : {};

/* destructor */
template <typename T>
mosquito_habitat<T>::~mosquito_habitat(){};

/* default move semantics */
template <typename T>
mosquito_habitat<T>::mosquito_habitat(mosquito_habitat&& rhs) = default;

template <typename T> mosquito_habitat<T>&
mosquito_habitat<T>::operator=(mosquito_habitat&& rhs) = default;



/* stochastic Euler-step */
template <>
inline void mosquito_habitat::euler_step<int>(const double dt){

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # EARLY-STAGE LARVAL INSTARS (EL)
  ######################################## */



  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

};


/* deterministic Euler-step */
template <>
inline void mosquito_habitat::euler_step<int>(const double dt){

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # EARLY-STAGE LARVAL INSTARS (EL)
  ######################################## */

  /* inbound oviposition to EL */
  double NV = SV + EV + IV;
  EL_new = (betaCom * NV * dt);

  /* instantaneous hazards for EL */
  double haz_EL_mort <- pars["muEL"]*(1 + ((EL+LL)/K));
  double haz_EL_2LL <- 1.0 / pars["durEL"];

  /* jump probabilities */
  EL_probs[0] = exp(-(haz_EL_mort + haz_EL_2LL)*dt);
  EL_probs[1] = (1 - EL_probs[0])*(haz_EL_mort / (haz_EL_mort + haz_EL_2LL)); /* death */
  EL_probs[2] = (1 - EL_probs[0])*(haz_EL_2LL / (haz_EL_mort + haz_EL_2LL)); /* to late-instar */

  /* jump sizes */
  EL_transitions[0] = EL * EL_probs[0];
  EL_transitions[1] = EL * EL_probs[1];
  EL_transitions[2] = EL * EL_probs[2];

  /* ########################################
  # LATE-STAGE LARVAL INSTARS (LL)
  ######################################## */

  /* instantaneous hazards for LL */
  double haz_LL_mort = pars["muLL"]*(1.0 + pars["gamma"]*((EL+LL)/K));
  double haz_LL_2PL = 1.0 / pars["durLL"];

  /* jump probabilities */
  LL_probs[0] = exp(-(haz_LL_mort + haz_LL_2PL)*dt);
  LL_probs[1] = (1 - LL_probs[0])*(haz_LL_mort / (haz_LL_mort + haz_LL_2PL)); /* death */
  LL_probs[2] = (1 - LL_probs[0])*(haz_LL_2PL / (haz_LL_mort + haz_LL_2PL)); /* to pupae */

  /* jump sizes */
  LL_transitions[0] = LL * LL_probs[0];
  LL_transitions[1] = LL * LL_probs[]1;
  LL_transitions[2] = LL * LL_probs[]2;

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

};




// template <typename T>
// class summary_base {
// public:
//   summary_base();               /* constructor */
//   virtual ~summary_base() = 0;  /* pure virtual destructor */
//
//   /*
//    * all summary classes must return a SEXP upon request (when the simulation is done)
//    * that will be returned to R
//    * (child classes inherit 'interface-only')
//    */
//   virtual SEXP output() = 0;
//
//   /* delete all copy semantics */
//   summary_base(const summary_base&) = delete;
//   summary_base& operator=(const summary_base&) = delete;
//
//   /* default move semantics */
//   summary_base(summary_base&&);
//   summary_base& operator=(summary_base&&);
// };
