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
void house::update_intervention(){
  if(IRS && village_ptr->get_tNow() > IRSoff){
    IRS = false;
  }
};

bool house::has_IRS(){
  return IRS;
}

void house::apply_IRS(){
  IRS = true;
  IRSoff = R::rexp(village_ptr->param_ptr->at("IRSduration"));
}

/* humans */
void house::distribute_EIR(){

  // DEBUG
  if(pi.size() != id.size() || id.size() != humans.size()) {
    Rcpp::Rcout << "error: size of id, pi, and humans vec is not all the same" << std::endl;
  }
  // DEBUG

  /* normalize pi */
  normalize_pi();

  /* only partition out bites if there's (attempted) bites to give */
  if(EIR > 0){
    /* if there's only one person, they get all the bites */
    if(humans.size() == 1){
      humans.front()->get_bitten(EIR);
    /* if there's multiple people, we have to partition bites accordingly */
    } else {

      /* rmultinom(int n, double* prob, int K, int* rN) */
      size_t K = pi.size();
      std::vector<int> EIR_2humans(K,0);
      rmultinom(EIR, pi.data(), K, EIR_2humans.data());

      /* sent the sampled bites to the human class */
      for(size_t k=0; k<K; k++){
        humans[k]->get_bitten(EIR_2humans[k]);
      }

    }
  /* no bites showed up here */
  } else {
    for(auto& h : humans){
      h->get_bitten(0);
    }
  }

};

void house::add_human(human_ptr h){
  id.emplace_back(h->get_id());
  pi.emplace_back(h->get_pi());
  humans.emplace_back(std::move(h));
  // // if this is a bug, pull out the ID first and pass to the std::pair obj
  // pi.emplace(h->get_id(),h->get_pi());
  // humans.emplace(h->get_id(),std::move(h));
};

void house::normalize_pi(){

  /* sum it */
  const double pi_sum = std::accumulate(pi.begin(), pi.end(), 0.0,
                                             [](const double previous, const double element)
                                              { return previous + element; });
  /* normalize it */
  std::for_each(pi.begin(), pi.end(), [pi_sum](double& element){
    element /= pi_sum;
  });
  /* bop it */

  // /* sum it */
  // const double pi_sum = std::accumulate(pi.begin(), pi.end(), 0.0,
  //                                            [](const double previous, const auto& element)
  //                                             { return previous + element.second; });
  // /* normalize it */
  // std::for_each(pi.begin(), pi.end(), [pi_sum](auto& element){
  //   element.second /= pi_sum;
  // });
  // /* bop it */
}
