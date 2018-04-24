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
 #  House Class Definition
*/

#ifndef RACD_HOUSE_
#define RACD_HOUSE_

#include <memory>
#include <utility>
#include <vector>

// #include "DEBUG.hpp"

/* typedefs and forward declarations */
class human;
typedef std::unique_ptr<human>          human_ptr;
typedef std::vector<human_ptr>          human_vector;

/* house */
class house {

public:

  /* constructor */
  house(const int& _houseID, const double& _psi, const double& _x, const double& _y);

  /* destructor */
  ~house();

  /* add humans */
  void                                      add_human(human_ptr h);

  /* accessors */
  int                                       get_houseID(){ return houseID; };
  double                                    get_psi(){ return psi; };
  human_vector&                             get_humans(){ return humans; };
  double                                    get_x(){ return x; };
  double                                    get_y(){ return y; };

private:

  int                                       houseID; /* ID */
  double                                    psi; /* relative risk */

  double                                    x; /* longitude */
  double                                    y; /* latitude */

  human_vector                              humans; /* people here */

};

#endif
