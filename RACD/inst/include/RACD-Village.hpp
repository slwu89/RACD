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
#include <string>

/* typedefs and forward declarations */
class house;
typedef std::unique_ptr<house> house_ptr;

class village {

public:

  /* constructor */
  village();

  /* destructor */
  ~village();

  /* demographics */
  void                                      births();


private:

  std::vector<house_ptr>                    houses;
  // std::vector<breedingSite_ptr>             breedingSites;

  std::string                               output; /* string giving valid output file to make */

};

#endif
