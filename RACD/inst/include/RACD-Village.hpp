/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  December 2017
 #
 #  Village Class Declaration
*/

/* ######################################################################
 # includes and foward declarations
###################################################################### */

#ifndef RACD_VILLAGE
#define RACD_VILLAGE

/* C++ includes */
#include <memory>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <math.h> /* sqrt */

/* Rcpp includes */
#include <Rcpp.h>

// #include "DEBUG.hpp"

/* alias and forward declarations */
class house;
using house_ptr = std::unique_ptr<house>;

using parameters = std::unordered_map<std::string,double>;

class mosquito_habitat;
using mosy_ptr = std::unique_ptr<mosquito_habitat>;


/* ######################################################################
 # class declaration
###################################################################### */

class village {

public:

  /* constructor & destructor */
  village(const Rcpp::List& theta, const Rcpp::IntegerVector& mosy);
  ~village();

  /* initialize objects */
  void initialize(const Rcpp::List& humansR, const Rcpp::List& housesR);

  /* Simulation Methods */

  /* one simulation run */
  // void                                      simulation(const int tMax);

  /* time */
  int                                       tNow;
  int                                       get_tNow(){return tNow;}

  /* daily simulation */
  void                                      one_day();

  /* demographics */
  void                                      births();
  void                                      deaths();

  /* utility classes */
  std::unique_ptr<parameters>               param_ptr;

  /* the landscape */
  mosy_ptr                                  mosquito;
  std::vector<house_ptr>                    houses;

  /* track output */
  void                                      track_human_state(Rcpp::IntegerMatrix& out);
  void                                      track_mosquito_state(Rcpp::IntegerMatrix& out);

private:
  /* keep track of all human IDs */
  int                                       max_humanID;

  /* for output */
  Rcpp::IntegerVector human_state;

};

#endif
