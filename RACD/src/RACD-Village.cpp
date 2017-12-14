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
 #  Village Class Implementation
*/

#include "RACD-Village.hpp"
#include "RACD-House.hpp"

/* constructor */
village::village(const Rcpp::List& humans, const Rcpp::List& houses){
  #ifdef DEBUG_HPP
  std::cout << "village being born at " << this << std::endl;
  #endif
};

/* destructor */
village::~village(){
  #ifdef DEBUG_HPP
  std::cout << "village being killed at " << this << std::endl;
  #endif
};

// /* Simulation Methods */
//
// /* daily simulation */
// void                                      one_day();
//
// /* demographics */
// void                                      births();
