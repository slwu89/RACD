#include "globals.hpp"

extern Rcpp::NumericMatrix output1;

// [[Rcpp::export]]
Rcpp::NumericMatrix test(){
  output1.fill(0);
  output1.at(0,0) += R::rpois(42);
  return output1;
}
