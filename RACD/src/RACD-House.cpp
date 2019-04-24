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
house::house(const int houseID_, village* const village_ptr_) : village_ptr(village_ptr_), houseID(houseID_),
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
  if(IRS && village_ptr->get_tNow() > IRSoff){
    IRS = false;
  }
  return IRS;
}

void house::apply_IRS(){
  IRS = true;
  IRSoff = R::rexp(village_ptr->param_ptr->at("IRSduration"));
}

/* humans */
void house::add_human(human_ptr h){
  // if this is a bug, pull out the ID first and pass to the std::pair obj
  pi.emplace(h->get_id(),h->get_pi());
  humans.emplace(h->get_id(),std::move(h));

  /* normalize pi */
  normalize_pi();
};

void house::normalize_pi(){
  /* sum it */
  const double pi_sum = std::accumulate(pi.begin(), pi.end(), 0.0,
                                             [](const double previous, const auto& element)
                                              { return previous + element.second; });
  /* normalize it */
  std::for_each(pi.begin(), pi.end(), [pi_sum](auto& element){
    element.second /= pi_sum;
  });
}
