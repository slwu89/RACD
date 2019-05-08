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
  double                                  psi;  // P(a bite going anywhere goes here)
  double                                  w;    // E[P(feed and survive)]
  double                                  y;    // E[P(feed)]
  double                                  z;    // E[P(repelled without feeding)]
  double                                  C;    // E[P(a feed here will result in a mosquito infection)]
  size_t                                  n;    // number of people here
  std::unordered_map<size_t,double>       pi;   // normalized PMF of who gets bitten
  double                                  EIR;  // the number of bites this house gets today
  bool                                    IRS;  // does my house have IRS

  std::list<human_ptr>                    humans; // who lives here

  // constructor & destructor
  house(const size_t id_, const double psi_, const double w_, const double y_, const double z_, const double C_,
        const size_t n_
  );
  ~house();
} house;

// pointer to house
using house_ptr = std::unique_ptr<house>;


/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh);

// normalize pi vector in a house
void normalize_pi(house_ptr& hh);

#endif
