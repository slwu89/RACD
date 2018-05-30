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
 #  PRNG Singleton
*/

#ifndef RACD_PRNG_
#define RACD_PRNG_

#include <random>
#include <iostream>

// #include "DEBUG.hpp"

/* threadsafe prng singleton */
class prng final {
public:
    /* utility methods */
    static prng&                           instance(); /* get instance */
    void                                   set_seed(const uint_least32_t& seed);

    /* continuous random variate sampling */
    double                                 get_runif();
    double                                 get_rexp(const double& rate);
    double                                 get_rlnorm(const double& meanlog, const double& sdlog);

    /* discrete random variate sampling */
    int                                    get_rpois(const double& lambda);
    int                                    get_rbinom(const int& n, const double& p);

    /* resample template type T x 'size' times */
    template<typename T>
    T                                      get_resample(const T& x, const int& size);

private:
  /* constructor & destructor */
  prng();
  ~prng();

  /* delete all copy & move semantics */
  prng(const prng&) = delete;
  prng& operator=(const prng&) = delete;
  prng(prng&&) = delete;
  prng& operator=(prng&&) = delete;

  std::mt19937                            rng;
  std::uniform_real_distribution<double>  runif;
};

#endif
