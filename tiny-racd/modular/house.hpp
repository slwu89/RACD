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

#ifndef HOUSE_HPP
#define HOUSE_HPP

#include <memory>
#include <list>
#include <unordered_map>
#include <algorithm>

#include <Rcpp.h>

// declare the human
struct human;
using human_ptr = std::unique_ptr<human>;


/* ################################################################################
#   the house
################################################################################ */

// a house
// all expectations are averaging out heterogeneity in the people at this house
// the reason pi is a hash table is because the humans are in a linked list
// so they will need to look up their pi by their ID number
typedef struct house {

  static int                              global_id;

  // data members
  int                                     id;
  double                                  psi;  // P(a bite going anywhere goes here)
  double                                  W;    // E[P(feed and survive)]
  double                                  Y;    // E[P(feed)]
  double                                  Z;    // E[P(repelled without feeding)]
  double                                  C;    // E[P(a feed here will result in a mosquito infection)]
  int                                     n;    // number of people here
  std::unordered_map<int,double>          pi;   // normalized PMF of who gets bitten
  double                                  EIR;  // the number of bites this house gets today

  // interventions
  bool                                    IRS;  // does my house have IRS
  int                                     IRS_time_off;
  int                                     cinc; // how many clinical cases happened today

  std::list<human_ptr>                    humans; // who lives here

  // constructor & destructor
  house();
  ~house();
} house;

int house::global_id = 0;

// pointer to house
using house_ptr = std::unique_ptr<house>;

// the houses
using house_vector = std::vector<house_ptr>;


/* ################################################################################
#   initialize houses and return to R
################################################################################ */

Rcpp::XPtr<house_vector> init_houses(
  const Rcpp::List& humans_param,
  const Rcpp::List& house_param
);


/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

// update global interface for mosquitos: THIS IS THE FUNCTION TO CALL (others are helpers)
// HUMAN -> MOSY
// returns bite_probs: vector of WW,ZZ,CC for mosquitoes
std::vector<double> update_biting(house_vector& houses);

// MOSY -> HUMAN (after running feeding_cycle)
// EIR is the vector from mosquitoes
void update_EIR(house_vector& houses, const std::vector<double>& EIR);

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh);

// update W (net probability of successfuly feeding)
void update_W(house_ptr& hh);

// update Z (net probability of being repelled w/out feeding)
void update_Z(house_ptr& hh);

// normalize pi vector in a house
void normalize_pi(house_ptr& hh);


/* ################################################################################
#   Interventions
################################################################################ */

// called in one_day_update
void update_interventions_house(house_ptr& hh, const int tnow);

// spray the house
void apply_IRS(house_ptr& hh, const int tnow);

// give ITNs to everyone in the house
void apply_ITN(house_ptr& hh);

// apply MDA; give drugs to everyone in the house indiscriminately
void apply_MDA(house_ptr& hh);

// apply RACD; test people via PCR, only give drugs to those who test positive
void apply_RACD_PCR(house_ptr& hh);

// apply RACD; test people via microscopy, only give drugs to those who test positive
void apply_RACD_Mic(house_ptr& hh);


/* ################################################################################
#   demographics
################################################################################ */

/* shamelessly "referenced" from Knuth, The Art of Computer Programming Vol 2, section 4.2.2 */
// the mean immunity among 18-22 year olds
double mean_ICA18_22(house_vector& houses);

// bring out yer dead!
void one_day_deaths(house_vector& houses);

// the respawn point
void one_day_births(house_vector& houses);


/* ################################################################################
#   daily simulation
################################################################################ */

// update the dynamics of a single dwelling
void one_day_update_house(house_ptr& hh, const int tnow);

// the update for all humans/dwellings
void one_day_update(house_ptr& houses, const int tnow);

// exposed to R
// runs one_day_update, one_day_births, one_day_deaths in that order
void one_day_step(SEXP houses, const int tnow);


/* ################################################################################
#   track output
################################################################################ */

Rcpp::DataFrame       track_state(SEXP houses, const int tnow);
Rcpp::DataFrame       track_immunity(SEXP houses, const int tnow);
Rcpp::DataFrame       track_transmission(SEXP houses, const int tnow);


#endif
