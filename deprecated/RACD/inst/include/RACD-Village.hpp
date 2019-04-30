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
#include <vector>
#include <iostream>
#include <math.h> /* sqrt */

#include <string>
#include <unordered_map>

/* Rcpp includes */
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// #include "DEBUG.hpp"

/* alias and forward declarations */
class house;
using house_ptr = std::unique_ptr<house>;

/* forward declare utility classes */
class prng;
class logger;
class parameters;


/* ######################################################################
 # class declaration
###################################################################### */

class village {

public:

  /* constructor & destructor */
  village(const uint_least32_t seed, const Rcpp::NumericVector& theta);
  ~village();

  /* initialize objects */
  void initialize(const Rcpp::List& humansR, const Rcpp::List& housesR);

  /* Simulation Methods */

  /* one simulation run */
  void                                      simulation(const int tMax, Rcpp::IntegerMatrix& state_out);

  /* daily simulation */
  void                                      one_day();

  /* demographics */
  void                                      births();
  void                                      deaths();

  /* utility classes */
  std::unique_ptr<prng>                     prng_ptr;
  // std::unique_ptr<logger>                   logger_ptr;
  std::unique_ptr<parameters>               param_ptr;

private:

  /* time */
  int                                       tNow;
  size_t                                    run_id;

  /* keep track of all human IDs */
  int                                       max_humanID;

  /* houses */
  std::vector<house_ptr>                    houses;
  // std::vector<breedingSite_ptr>             breedingSites;

};

#endif
