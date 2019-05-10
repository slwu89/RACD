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

// // output: states
// extern std::vector<size_t> state_S;
// extern std::vector<size_t> state_E;
// extern std::vector<size_t> state_T;
// extern std::vector<size_t> state_D;
// extern std::vector<size_t> state_A;
// extern std::vector<size_t> state_U;
// extern std::vector<size_t> state_P;
//
// // output: pop sizes
// extern std::vector<size_t> num_All;
// extern std::vector<size_t> num_2_10;
// extern std::vector<size_t> num_0_5;
// extern std::vector<size_t> num_5_10;
// extern std::vector<size_t> num_10_15;
// extern std::vector<size_t> num_15Plus;
//
// // output: clinical incidence
// extern std::vector<size_t> cinc_All;
// extern std::vector<size_t> cinc_2_10;
// extern std::vector<size_t> cinc_0_5;
// extern std::vector<size_t> cinc_5_10;
// extern std::vector<size_t> cinc_10_15;
// extern std::vector<size_t> cinc_15Plus;
//
// // output for mosquitos
// extern std::vector<size_t> mosy_S;
// extern std::vector<size_t> mosy_E;
// extern std::vector<size_t> mosy_I;
//
// // parameters
// extern std::unordered_map<std::string,double> parameters;
//
// // current simulation time
// extern size_t tnow;
//
// // id for people
// extern size_t global_hid;
//
// // transmission: mosy -> human
//
// // psi (biting weights) and EIR for each house
// extern std::vector<double>   psi;
// extern std::vector<double>   EIR;
//
// // transmission: human -> mosy
//
// extern double      CC; // P(bite will cause infection in mosquito) --- expectation of this prob over all landscape/individual heterogeneities
// extern double      WW; // avg probability to bite and survive
// extern double      ZZ; // avg probability to bite
//
// extern void reset_states(const size_t tmax);
// extern void reset_pop(const size_t tmax);
// extern void reset_cinc(const size_t tmax);
// extern void reset_mosy(const size_t tmax);
// extern void reset_globals(const size_t tmax);

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

  // constructor & destructor
  house(const size_t id_, const double psi_, const double W_, const double Y_, const double Z_, const double C_,
        const size_t n_
  );
  ~house();
} house;

// pointer to house
using house_ptr = std::unique_ptr<house>;

// the houses
using house_vector = std::vector<house_ptr>;


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

// update global interface for mosquitos
void update_biting(house_vector& houses);


/* ################################################################################
#   Interventions
################################################################################ */

void update_interventions(house_ptr& hh);

void apply_IRS(house_ptr& hh);

#endif
