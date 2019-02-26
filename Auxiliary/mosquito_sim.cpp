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

/* template for int and double */
template <typename T>
class mosquito_habitat {
public:
  mosquito_habitat();
  ~mosquito_habitat();

  /* simulation */
  void euler_step();

  

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
