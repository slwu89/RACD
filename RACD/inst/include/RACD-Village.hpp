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
 #  Village Class Definition
*/

#ifndef _RACD_VILLAGE_
#define _RACD_VILLAGE_

#include <Rcpp.h>
#include <memory>
#include <vector>

/* typedefs and forward declarations */
class house;
typedef std::unique_ptr<house> house_ptr;

class village {

public:

  /* constructor */
  village(const Rcpp::List& humans, const Rcpp::List& houses);

  /* destructor */
  ~village();

  /* Simulation Methods */

  /* daily simulation */
  void                                      one_day();

  /* demographics */
  void                                      births();


private:

  std::vector<house_ptr>                    houses;
  // std::vector<breedingSite_ptr>             breedingSites;


};

#endif
