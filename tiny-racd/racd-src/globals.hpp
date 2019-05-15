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
#include <memory>

// // all globals as struct
// typedef struct globals {
//
//   // output for humans
//
//   // output: states
//   std::vector<size_t> state_S;
//   std::vector<size_t> state_E;
//   std::vector<size_t> state_T;
//   std::vector<size_t> state_D;
//   std::vector<size_t> state_A;
//   std::vector<size_t> state_U;
//   std::vector<size_t> state_P;
//
//   // output: pop sizes
//   std::vector<size_t> num_All;
//   std::vector<size_t> num_2_10;
//   std::vector<size_t> num_0_5;
//   std::vector<size_t> num_5_10;
//   std::vector<size_t> num_10_15;
//   std::vector<size_t> num_15Plus;
//
//   // output: clinical incidence
//   std::vector<size_t> cinc_All;
//   std::vector<size_t> cinc_2_10;
//   std::vector<size_t> cinc_0_5;
//   std::vector<size_t> cinc_5_10;
//   std::vector<size_t> cinc_10_15;
//   std::vector<size_t> cinc_15Plus;
//
//   // output for mosquitos
//   std::vector<size_t> mosy_S;
//   std::vector<size_t> mosy_E;
//   std::vector<size_t> mosy_I;
//
//   // parameters
//   std::unordered_map<std::string,double> parameters;
//
//   // state variables
//
//   // current simulation time
//   size_t tnow;
//
//   // id for people
//   size_t global_hid;
//
//   // transmission: mosy -> human
//
//   // psi (biting weights) and EIR for each house
//   std::vector<double>   psi;
//   std::vector<double>   EIR;
//
//   // transmission: human -> mosy
//
//   double      CC; // P(bite will cause infection in mosquito) --- expectation of this prob over all landscape/individual heterogeneities
//   double      WW; // avg probability to bite and survive
//   double      ZZ; // avg probability to bite
//
//   // constructor & destructor
//   globals(const size_t tmax, const size_t nhouse);
//
//   ~globals();
// } globals;
//
// // constructor
// globals::globals(const size_t tmax, const size_t nhouse) :
//   state_S(tmax,0), state_E(tmax,0), state_T(tmax,0), state_D(tmax,0), state_A(tmax,0), state_U(tmax,0), state_P(tmax,0),
//   num_All(tmax,0), num_2_10(tmax,0), num_0_5(tmax,0), num_5_10(tmax,0), num_10_15(tmax,0), num_15Plus(tmax,0),
//   cinc_All(tmax,0), cinc_2_10(tmax,0), cinc_0_5(tmax,0), cinc_5_10(tmax,0), cinc_10_15(tmax,0), cinc_15Plus(tmax,0),
//   mosy_S(tmax,0), mosy_E(tmax,0), mosy_I(tmax,0),
//   parameters(), tnow(0), global_hid(0),
//   psi(nhouse,0.), EIR(nhouse,0.),
//   CC(0.), WW(0.), ZZ(0.)
// {};
//
// // destructor
// globals::~globals(){};
//
// // alias the ptr
// using globals_ptr = std::unique_ptr<globals>;


// output for humans

// output: states
extern std::vector<size_t> state_S;
extern std::vector<size_t> state_E;
extern std::vector<size_t> state_T;
extern std::vector<size_t> state_D;
extern std::vector<size_t> state_A;
extern std::vector<size_t> state_U;
extern std::vector<size_t> state_P;

// output: pop sizes
extern std::vector<size_t> num_All;
extern std::vector<size_t> num_2_10;
extern std::vector<size_t> num_0_5;
extern std::vector<size_t> num_5_10;
extern std::vector<size_t> num_10_15;
extern std::vector<size_t> num_15Plus;

// output: clinical incidence
extern std::vector<size_t> cinc_All;
extern std::vector<size_t> cinc_2_10;
extern std::vector<size_t> cinc_0_5;
extern std::vector<size_t> cinc_5_10;
extern std::vector<size_t> cinc_10_15;
extern std::vector<size_t> cinc_15Plus;

// output for mosquitos
extern std::vector<size_t> mosy_S;
extern std::vector<size_t> mosy_E;
extern std::vector<size_t> mosy_I;

// other output
extern std::vector<size_t> time_out;

// parameters
extern std::unordered_map<std::string,double> parameters;

// state variables

// current simulation time
extern size_t tnow;

// id for people
extern size_t global_hid;

// transmission: mosy -> human

// psi (biting weights) and EIR for each house
extern std::vector<double>   GLOBAL_PSI;
extern std::vector<double>   EIR;

extern size_t NHOUSE;

// transmission: human -> mosy

extern double      CC; // P(bite will cause infection in mosquito) --- expectation of this prob over all landscape/individual heterogeneities
extern double      WW; // avg probability to bite and survive
extern double      ZZ; // avg probability to bite


// function to reset globals between calls from R
inline void reset_states(const size_t tmax){
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

  time_out.clear();
  time_out.resize(tmax,0);
}

inline void reset_pop(const size_t tmax){
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

inline void reset_cinc(const size_t tmax){
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

inline void reset_mosy(const size_t tmax){
  mosy_S.clear();
  mosy_E.clear();
  mosy_I.clear();
  mosy_S.resize(tmax,0);
  mosy_E.resize(tmax,0);
  mosy_I.resize(tmax,0);
}

// reset all of the things
inline void reset_globals(const size_t tmax, const size_t nhouse){

  parameters.clear();

  GLOBAL_PSI.clear();
  EIR.clear();
  GLOBAL_PSI.resize(nhouse,0);
  EIR.resize(nhouse,0);

  NHOUSE = nhouse;

  tnow = 0;
  global_hid = 0;

  CC = 0.;
  WW = 0.;
  ZZ = 0.;

  reset_states(tmax);
  reset_pop(tmax);
  reset_cinc(tmax);
  reset_mosy(tmax);

};

#endif
