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
#include <list>

/* typedefs and forward declarations */
class human;
typedef std::unique_ptr<human> human_ptr;

/* house */
class house {

public:

  house();
  ~house();


private:

  int                                       houseID; /* ID */

  double                                    x; /* longitude */
  double                                    y; /* latitude */

  std::list<human_ptr>                      humans; /* people here */

};

#endif
