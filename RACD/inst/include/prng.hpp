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

#include "debug.hpp"

class prng {
public:
    prng(const uint_least32_t &seed) : rng(seed) {
      runif = std::uniform_real_distribution<double>(0,1);
      #ifdef DEBUG_HPP
      std::cout << "prng being born at " << this << std::endl;
      #endif
    };
    ~prng(){
      #ifdef DEBUG_HPP
      std::cout << "prng being killed at " << this << std::endl;
      #endif
    };

    double get_runif(){
        return runif(rng);
    };

private:
    std::mt19937                            rng;
    std::uniform_real_distribution<double>  runif;

};

#endif
