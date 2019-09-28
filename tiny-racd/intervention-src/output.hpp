/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu
 #  slwu89@berkeley.edu
 #  September 2019
 #
 #  stuff the whole program needs to see
*/

#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <vector>

#include <Rcpp.h>


class output final {
public:
  /* utility methods */
  static output&                           instance(); /* get instance */

private:

  /* constructor & destructor */
  output();
  ~output();

  /* delete all copy & move semantics */
  output(const output&) = delete;
  output& operator=(const output&) = delete;
  output(output&&) = delete;
  output& operator=(output&&) = delete;

  /* output */
  Rcpp::List              state_age; // states X age
  Rcpp::IntegerMatrix     cinc_age; // clinical incidence X age
  Rcpp::IntegerMatrix     mosquito; // SEI

};

#endif
