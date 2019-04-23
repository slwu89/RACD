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
#include "RACD-Village.hpp"


/* constructor */
house::house(const int houseID_, const double psi_, village* const village_ptr_) : village_ptr(village_ptr_), houseID(houseID_), psi(psi_),
             EIR(0.0), IRS(false), IRSoff(2E16)
{
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

/* interventions */
bool house::has_IRS(){
  if(!IRS){
    return false;
  } else {
    if(village_ptr->get_tNow() > IRSoff){
      IRS = false;
      return false;
    } else {
      return true;
    }
  }
}

void house::apply_IRS(){
  IRS = true;
  IRSoff = R::rexp(village_ptr->param_ptr->at("IRSduration"));
}

/* add humans */
void house::add_human(human_ptr h){
  humans.push_back(std::move(h));
};
