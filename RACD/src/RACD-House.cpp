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
 #  House Class Implementation
*/

#include "RACD-House.hpp"
#include "RACD-Human.hpp"


/* constructor */
house::house(const int& _houseID, const double& _psi) : houseID(_houseID), psi(_psi) {
  #ifdef DEBUG_HPP
  std::cout << "house " << houseID << " being born at " << this << std::endl;
  #endif
};

/* destructor */
house::~house(){
  #ifdef DEBUG_HPP
  std::cout << "house " << houseID << " being killed at " << this << std::endl;
  #endif
};
