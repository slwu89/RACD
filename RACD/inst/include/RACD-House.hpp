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
  house();

  /* destructor */
  ~house();

  /* accessors */
  int                                       get_houseID(){ return houseID; };
  human_vector&                             get_humans(){ return humans; };


private:

  int                                       houseID; /* ID */

  double                                    x; /* longitude */
  double                                    y; /* latitude */

  human_vector                              humans; /* people here */

};

#endif
