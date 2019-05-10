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
void reset_globals(){
  Rcpp::Rcout << "write me\n";
};


#endif
