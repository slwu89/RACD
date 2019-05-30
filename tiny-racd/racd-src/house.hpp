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

#ifndef house_hpp
#define house_hpp

#include <Rcpp.h>

#include <memory>
#include <list>
#include <unordered_map>
#include <algorithm>

// gloabl online stats
class RunningStat;

// declare the human
struct human;
using human_ptr = std::unique_ptr<human>;

// a house
// all expectations are averaging out heterogeneity in the people at this house
// the reason pi is a hash table is because the humans are in a linked list
// so they will need to look up their pi by their ID number
typedef struct house {

  // data members
  size_t                                  id;
  // double                                  psi;  // P(a bite going anywhere goes here)
  double                                  W;    // E[P(feed and survive)]
  double                                  Y;    // E[P(feed)]
  double                                  Z;    // E[P(repelled without feeding)]
  double                                  C;    // E[P(a feed here will result in a mosquito infection)]
  size_t                                  n;    // number of people here
  std::unordered_map<size_t,double>       pi;   // normalized PMF of who gets bitten
  double                                  EIR;  // the number of bites this house gets today

  // interventions
  bool                                    IRS;  // does my house have IRS
  size_t                                  IRS_time_off;

  std::list<human_ptr>                    humans; // who lives here

  // stats
  RunningStat* const                      global_stat;

  // constructor & destructor
  house(const size_t id_, RunningStat* global_stat_);

  ~house();
} house;

// pointer to house
using house_ptr = std::unique_ptr<house>;

// the houses
using house_vector = std::vector<house_ptr>;


/* ################################################################################
#   tracking output
################################################################################ */

void track_state(const house_vector& houses);

void track_age(const house_vector& houses);


/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh);

// update W (net probability of successfuly feeding)
void update_W(house_ptr& hh);

// update Z (net probability of being repelled w/out feeding)
void update_Z(house_ptr& hh);

// normalize pi vector in a house
void normalize_pi(house_ptr& hh);

// update global interface for mosquitos: THIS IS THE FUNCTION TO CALL (others are helpers)
// HUMAN -> MOSY
void update_biting(house_vector& houses);

// MOSY -> HUMAN (after running feeding_cycle)
void update_EIR(house_vector& houses);


/* ################################################################################
#   Interventions
################################################################################ */

void update_interventions_house(house_ptr& hh);

void apply_IRS(house_ptr& hh);


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
void one_day_update_house(house_ptr& hh);

// the update for all humans/dwellings
void one_day_update(house_vector& houses);

#endif
