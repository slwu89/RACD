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
 #  stuff the whole program needs to see
*/

#ifndef globals_hpp
#define globals_hpp

#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>

// output for humans

// output: states
std::vector<size_t> state_S(1,0);
std::vector<size_t> state_E(1,0);
std::vector<size_t> state_T(1,0);
std::vector<size_t> state_D(1,0);
std::vector<size_t> state_A(1,0);
std::vector<size_t> state_U(1,0);
std::vector<size_t> state_P(1,0);

// output: pop sizes
std::vector<size_t> num_All(1,0);
std::vector<size_t> num_2_10(1,0);
std::vector<size_t> num_0_5(1,0);
std::vector<size_t> num_5_10(1,0);
std::vector<size_t> num_10_15(1,0);
std::vector<size_t> num_15Plus(1,0);

// output: clinical incidence
std::vector<size_t> cinc_All(1,0);
std::vector<size_t> cinc_2_10(1,0);
std::vector<size_t> cinc_0_5(1,0);
std::vector<size_t> cinc_5_10(1,0);
std::vector<size_t> cinc_10_15(1,0);
std::vector<size_t> cinc_15Plus(1,0);

// output for mosquitos
std::vector<size_t> mosy_S(1,0);
std::vector<size_t> mosy_E(1,0);
std::vector<size_t> mosy_I(1,0);

// parameters
std::unordered_map<std::string,double> parameters;

// current simulation time
size_t tnow = 0;

// transmission: mosy -> human

// psi (biting weights) and EIR for each house
std::vector<double>   psi;
std::vector<double>   EIR;

// transmission: human -> mosy

double      CC; // P(bite will cause infection in mosquito) --- expectation of this prob over all landscape/individual heterogeneities
double      WW; // avg probability to bite and survive
double      ZZ; // avg probability to bite

// function to reset globals between calls from R
void reset_states(const size_t tmax){
  state_S.clear();
  state_E.clear();
  state_T.clear();
  state_D.clear();
  state_A.clear();
  state_U.clear();
  state_P.clear();
  state_S.resize(tmax,0);
  state_E.resize(tmax,0);
  state_T.resize(tmax,0);
  state_D.resize(tmax,0);
  state_A.resize(tmax,0);
  state_U.resize(tmax,0);
  state_P.resize(tmax,0);
}

void reset_pop(const size_t tmax){
  num_All.clear();
  num_2_10.clear();
  num_0_5.clear();
  num_5_10.clear();
  num_10_15.clear();
  num_15Plus.clear();
  num_All.resize(tmax,0);
  num_2_10.resize(tmax,0);
  num_0_5.resize(tmax,0);
  num_5_10.resize(tmax,0);
  num_10_15.resize(tmax,0);
  num_15Plus.resize(tmax,0);
};

void reset_cinc(const size_t tmax){
  cinc_All.clear();
  cinc_2_10.clear();
  cinc_0_5.clear();
  cinc_5_10.clear();
  cinc_10_15.clear();
  cinc_15Plus.clear();
  cinc_All.resize(tmax,0);
  cinc_2_10.resize(tmax,0);
  cinc_0_5.resize(tmax,0);
  cinc_5_10.resize(tmax,0);
  cinc_10_15.resize(tmax,0);
  cinc_15Plus.resize(tmax,0);
};

void reset_mosy(const size_t tmax){
  mosy_S.clear();
  mosy_E.clear();
  mosy_I.clear();
  mosy_S.resize(tmax,0);
  mosy_E.resize(tmax,0);
  mosy_I.resize(tmax,0);
}

// reset all of the things
void reset_globals(const size_t tmax){

  parameters.clear();
  psi.clear();
  EIR.clear();

  tnow = 0;

  CC = 0.;
  WW = 0.;
  ZZ = 0.;

  reset_states(tmax);
  reset_pop(tmax);
  reset_cinc(tmax);
  reset_mosy(tmax);

};

#endif
