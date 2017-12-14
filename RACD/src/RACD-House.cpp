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
house::house(const int& _houseID, const double& _psi, const double& _x, const double& _y) : houseID(_houseID), psi(_psi), x(_x), y(_y) {
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
