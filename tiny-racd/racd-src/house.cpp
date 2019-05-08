/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  April 2019
 #
 #  the house
*/

#include "house.hpp"

// other headers we need
#include "human.hpp"
#include "globals.hpp"


/* ################################################################################
#   constructor & destructor
################################################################################ */

house::house(const size_t id_, const double psi_, const double w_, const double y_, const double z_, const double C_,
             const size_t n_
) :
  id(id_), psi(psi_), w(w_), y(y_), z(z_), C(C_), n(n_), EIR(0), IRS(false)
{};

house::~house(){};


/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh){

  hh->C = 0.;

  for(auto& h : hh->humans){
    hh->C += hh->pi.at(h->id) * get_w(h) * h->c;
  }

};

// normalize pi vector in a house
void normalize_pi(house_ptr& hh){

  const double pi_sum = std::accumulate(hh->pi.begin(), hh->pi.end(), 0.0,
    [](const double prev, const auto& element){
      return prev + element.second;
    }
  );

  std::for_each(hh->pi.begin(), hh->pi.end(), [pi_sum](auto& element){
    element.second /= pi_sum;
  });


};
