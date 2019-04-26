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
 #  House Class Declaration
*/

/* ######################################################################
 # includes and foward declarations
###################################################################### */

#ifndef RACD_HOUSE
#define RACD_HOUSE

#include <memory>
#include <utility>
#include <algorithm>
#include <numeric> // for accumulate
#include <vector>

#include <unordered_map>

#include <Rcpp.h>
#include <Rmath.h> // for rmultinom

// #include "DEBUG.hpp"

/* alias and forward declarations */
class human;
using human_ptr       = std::unique_ptr<human>;

using human_vector    = std::vector<human_ptr>;
using human_pi        = std::vector<double>;
using human_id        = std::vector<int>;

class village;

/* for output */
static const std::unordered_map<std::string,size_t> state2col = {
  {"S",0},
  {"T",1},
  {"D",2},
  {"A",3},
  {"U",4},
  {"P",5}
};


/* ######################################################################
 # class declarations
###################################################################### */

class house {

public:

  /* constructor */
  house(const int houseID_, village* const village_ptr_);

  /* destructor */
  ~house();

  /* humans */
  void                                      add_human(human_ptr h);
  void                                      normalize_pi();

  /* accessors */
  int                                       houseID; /* ID */

  /* humans */
  human_vector                              humans; /* people here */
  human_pi                                  pi; /* biting weight on humans */
  human_id                                  id; /* biting weight on humans */

  /* biting */
  void                                      distribute_EIR();
  int                                       get_EIR(){return EIR;};
  void                                      set_EIR(const int e){ EIR = e;};

  /* interventions */
  void                                      update_intervention();

  bool                                      has_IRS();
  void                                      apply_IRS();

  /* track output */
  Rcpp::IntegerVector                       output_states();

  village* const                            village_ptr; /* raw pointer ok because village lifespan > house lifespan */

private:

  int                                       EIR; /* the total # of infectious bites arriving at this house, today */

  bool                                      IRS; /* does this house have IRS */
  double                                    IRSoff;

  // /* for output */
  // Rcpp::IntegerVector                       states;
};

#endif
