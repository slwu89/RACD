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
 #  Anomalous Materials Lab
*/

#include <Rcpp.h>
#include "parameters.hpp"

#include <iostream>

//' @export
// [[Rcpp::export]]
void test_parameters(){
  std::cout << "testing singleton parameters" << std::endl;
}