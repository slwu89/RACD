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

#ifndef _RACD_HOUSE_
#define _RACD_HOUSE_

#include <Rcpp.h>
#include <memory>
#include <vector>

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
