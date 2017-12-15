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
 #  PRNG Class Definition & Implementation
*/

#ifndef _RACD_PRNG_
#define _RACD_PRNG_

#include <random>
#include <iostream>

#include "DEBUG.hpp"

class prng {
public:

    /* constructor */
    prng(const uint_least32_t &seed) : rng(seed) {
      runif = std::uniform_real_distribution<double>(0,1);
      #ifdef DEBUG_HPP
      std::cout << "prng being born at " << this << std::endl;
      #endif
    };

    /* destructor */
    ~prng(){
      #ifdef DEBUG_HPP
      std::cout << "prng being killed at " << this << std::endl;
      #endif
    };

    /* runif(0,1) */
    double get_runif(){
      return runif(rng);
    };

    /* rbinom(n,p) */
    int get_rbinom(const int& n, const double& p){
      std::binomial_distribution<int>rbinom(n,p);
      return rbinom(rng);
    };

private:
    std::mt19937                            rng;
    std::uniform_real_distribution<double>  runif;

};

#endif
