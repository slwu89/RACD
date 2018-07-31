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
house::house(const int houseID_, const double psi_, const double x_, const double y_, village* village_ptr_) : houseID(houseID_), psi(psi_), x(x_), y(y_), village_ptr(village_ptr_) {
  #ifdef DEBUG_RACD
  std::cout << "house " << houseID << " being born at " << this << std::endl;
  #endif
};

/* destructor */
house::~house(){
  #ifdef DEBUG_RACD
  std::cout << "house " << houseID << " being killed at " << this << std::endl;
  #endif
};

/* add humans */
void house::add_human(human_ptr h){
  humans.push_back(std::move(h));
};
