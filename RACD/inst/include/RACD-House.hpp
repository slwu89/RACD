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
#include <unordered_map>
#include <algorithm>
#include <numeric> // for accumulate
#include <vector>

#include <Rmath.h> // for rmultinom

// #include "DEBUG.hpp"

/* alias and forward declarations */
class human;
using human_ptr       = std::unique_ptr<human>;

// using human_table     = std::unordered_map<int,human_ptr>;
// using human_pi        = std::unordered_map<int,double>;

using human_vector    = std::vector<human_ptr>;
using human_pi        = std::vector<double>;
using human_id        = std::vector<int>;

class village;


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
  int                                       get_houseID(){ return houseID; };


  // human_table&                              get_humans(){ return humans; };
  // double                                    get_pi(const int id){ return pi.at(id); };

  /* the humans */
  human_vector                                  humans; /* people here */
  human_pi                                     pi; /* biting weight on humans */
  human_id                                  id; /* biting weight on humans */

  /* biting */
  void                                      distribute_EIR();
  int                                       get_EIR(){return EIR;};
  void                                      set_EIR(const int e){ EIR = e;};

  /* interventions */
  void                                      update_intervention();

  bool                                      has_IRS();
  void                                      apply_IRS();

  village* const                            village_ptr; /* raw pointer ok because village lifespan > house lifespan */

private:

  int                                       houseID; /* ID */
  int                                       EIR; /* the total # of infectious bites arriving at this house, today */

  bool                                      IRS; /* does this house have IRS */
  double                                    IRSoff;
};

#endif
