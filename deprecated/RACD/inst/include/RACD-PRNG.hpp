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
 #  PRNG Class Declaration
*/

/* ######################################################################
 # includes and foward declarations
###################################################################### */

#ifndef RACD_PRNG
#define RACD_PRNG

/* C++ includes */
#include <random>
#include <iostream>

// #include "DEBUG.hpp"

/* ######################################################################
 # class declaration
###################################################################### */

class prng {
public:

  /* constructor & destructor */
  prng(const uint_least32_t seed);
  ~prng();

  /* delete all copy & move semantics */
  prng(const prng&) = delete;
  prng& operator=(const prng&) = delete;
  prng(prng&&) = delete;
  prng& operator=(prng&&) = delete;

  /* continuous random variate sampling */
  double                                 get_runif();
  double                                 get_rexp(const double& rate);
  double                                 get_rlnorm(const double& meanlog, const double& sdlog);
  double                                 get_beta_1_b(const double b){return 1.0 - pow(runif(rng), 1.0/b);};

  /* discrete random variate sampling */
  int                                    get_rpois(const double& lambda);
  int                                    get_rbinom(const int& n, const double& p);

  /* multivariate discrete sampling */
  void                                  get_rmultinom(int size, const std::vector<double>& prob, std::vector<int>& out, double switchover = 1.0);

  /* resample template type T x 'size' times */
  template<typename T>
  T                                      get_resample(const T& x, const int& size);

private:
  std::mt19937                            rng;
  std::uniform_real_distribution<double>  runif;
};

#endif
