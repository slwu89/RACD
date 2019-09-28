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

#include "globals.hpp"

// definitions of the things we declared extern in the header

// // age-stratified output
// Rcpp::IntegerMatrix state_age(7,6);
// Rcpp::rownames(state_age) = Rcpp::CharacterVector::create("S","E","T","D","A","U","P");
// Rcpp::colnames(state_age) = Rcpp::CharacterVector::create("all","2-10","0-5","5-10","10-15","15+");

// output: states
std::vector<size_t> state_S;
std::vector<size_t> state_E;
std::vector<size_t> state_T;
std::vector<size_t> state_D;
std::vector<size_t> state_A;
std::vector<size_t> state_U;
std::vector<size_t> state_P;

// output: pop sizes
std::vector<size_t> num_All;
std::vector<size_t> num_2_10;
std::vector<size_t> num_0_5;
std::vector<size_t> num_5_10;
std::vector<size_t> num_10_15;
std::vector<size_t> num_15Plus;

// output: clinical incidence
std::vector<size_t> cinc_All;
std::vector<size_t> cinc_2_10;
std::vector<size_t> cinc_0_5;
std::vector<size_t> cinc_5_10;
std::vector<size_t> cinc_10_15;
std::vector<size_t> cinc_15Plus;

// output for mosquitos
std::vector<size_t> mosy_S;
std::vector<size_t> mosy_E;
std::vector<size_t> mosy_I;

// transmision output
std::vector<double> lambda_v;
std::vector<double> eir_mean;
std::vector<double> eir_var;
// std::vector<double> lambda_h_mean;
// std::vector<double> lambda_h_var;
std::vector<double> b_mean;
std::vector<double> b_var;

// other output
std::vector<size_t> time_out;

// parameters
std::unordered_map<std::string,double> parameters;

// current simulation time
size_t tnow;

// id for people
size_t global_hid;

// transmission: mosy -> human

// psi (biting weights) and EIR for each house
std::vector<double>   psi;
std::vector<double>   EIR;

size_t NHOUSE;

// transmission: human -> mosy

double      CC; // P(bite will cause infection in mosquito) --- expectation of this prob over all landscape/individual heterogeneities
double      WW; // avg probability to bite and survive
double      ZZ; // avg probability to bite
