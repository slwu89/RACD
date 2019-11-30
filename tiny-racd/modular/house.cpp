/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  November 2019
 #
 #  Mosquitoes
*/

#include "house.hpp"
#include "human.hpp"

// [[Rcpp::plugins(cpp14)]]


/* ################################################################################
#   house ctor/dtor
################################################################################ */

house::house() :
  id(global_id++), W(0.), Y(0.), Z(0.), C(0.), n(0), EIR(0), IRS(false), IRS_time_off(0.), cinc(0)
{
  Rcpp::Rcout << "house id: " << id << ", ctor called at " << this << "\n";
};

house::~house(){
  Rcpp::Rcout << "house id: " << id << ", dtor called at " << this << "\n";
};
