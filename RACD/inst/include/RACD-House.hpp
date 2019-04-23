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
#include <vector>

// #include "DEBUG.hpp"

/* alias and forward declarations */
class human;
using human_ptr       = std::unique_ptr<human>;
using human_vector    = std::vector<human_ptr>;

class village;


/* ######################################################################
 # class declarations
###################################################################### */

class house {

public:

  /* constructor */
  house(const int houseID_, const double psi_, village* const village_ptr_);

  /* destructor */
  ~house();

  /* add humans */
  void                                      add_human(human_ptr h);

  /* accessors */
  int                                       get_houseID(){ return houseID; };
  double                                    get_psi(){ return psi; };
  human_vector&                             get_humans(){ return humans; };

  double                                    get_EIR(){return EIR;};
  void                                      set_EIR(const double e){ EIR = e;};

  bool                                      has_IRS();
  void                                      apply_IRS();

  village* const                            village_ptr; /* raw pointer ok because village lifespan > house lifespan */

private:

  int                                       houseID; /* ID */
  double                                    psi; /* relative risk */
  double                                    EIR; /* the total EIR at this house */

  bool                                      IRS; /* does this house have IRS */
  double                                    IRSoff;

  human_vector                              humans; /* people here */
  std::vector<double>                       pi; /* biting weight on humans */
};

#endif
