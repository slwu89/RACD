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

house::house(const size_t id_, const double psi_, const double W_, const double Y_, const double Z_, const double C_,
             const size_t n_
) :
  id(id_), psi(psi_), W(W_), Y(Y_), Z(Z_), C(C_), n(n_), EIR(0), IRS(false), IRS_time_off(0.)
{};

house::~house(){};


/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

// for each house

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh){

  hh->C = 0.;

  for(auto& h : hh->humans){
    hh->C += hh->pi.at(h->id) * get_w(h) * h->c;
  }

};

// update W (net probability of successfuly feeding)
void update_W(house_ptr& hh){

  hh->W = 0.;

  for(auto& h : hh->humans){
    hh->W += hh->pi.at(h->id) * get_w(h);
  }

};

// update Z (net probability of being repelled w/out feeding)
void update_Z(house_ptr& hh){

  hh->Z = 0.;

  for(auto& h : hh->humans){
    hh->Z += hh->pi.at(h->id) * get_z(h);
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


// for the global landscape

// this function updates (in order):
// pi for each house
// C for each house
// global CC
// global WW
// global ZZ
void update_biting(house_vector& houses){

  CC = 0.;
  WW = 0.;
  ZZ = 0.;
  for(auto& hh : houses){

    /* normalize biting weights (pi) */
    normalize_pi(hh);

    /* C for each house */
    update_C(hh);

    /* W for each house */
    update_W(hh);

    /* Z for each house */
    update_Z(hh);

    double psi_h = psi.at(hh->id);

    /* global CC is weighted average of household level C's */
    CC += psi_h * hh->C;

    /* global W is a weighted average of household level W's */
    WW += psi_h * hh->W;

    /* global Z is a weighted average of household level Z's */
    WW += psi_h * hh->Z;
  }

};


/* ################################################################################
#   Interventions
################################################################################ */

void update_interventions(house_ptr& hh){
  if(hh->IRS && tnow > hh->IRS_time_off){
    hh->IRS = false;
  }
};

void apply_IRS(house_ptr& hh){

  double IRS_decay = parameters.at("IRS_decay");

  hh->IRS = true;
  hh->IRS_time_off = R::rgeom(IRS_decay);

};
