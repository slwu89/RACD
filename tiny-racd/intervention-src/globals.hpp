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

#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include <vector>
#include <string>
#include <array>
#include <unordered_map>

#include <Rcpp.h>

class RunningStat;
using stat_map = std::unordered_map<std::string,RunningStat>;

using pmap = std::unordered_map<std::string,double>;
using dvec = std::vector<double>;


class globals final {
public:

  /* utility methods */
  static globals&         instance(); /* get instance */

  // set up size of output objects
  void                    set_output(const size_t tmax_);

  // call this when we're ready to return everything to R
  Rcpp::List              get_output();

  // set up the paramters hash table
  void                    set_parameters(const Rcpp::NumericVector& theta);

  // push the running statistics and clear before the next day
  void                    init_stats();
  stat_map&               get_stats(){return running_stats;}
  void                    push_stats();

  // called by mosquito::feeding_cycle
  void                    update_EIR(const double bites);

  // do this at the very bottom of the loop (just before the clock ticks to tomorrow)
  void                    iterate();

  // accessors
  void                    set_NHOUSE(const size_t NHOUSE_){NHOUSE = NHOUSE_;}
  size_t                  get_NHOUSE(){return NHOUSE;}

  pmap&                   get_pmap(){return parameters;}

  size_t                  get_tnow(){return tnow;}
  void                    inc_tnow(){tnow++;}
  bool                    tcheck(){return tnow < tmax;}

  size_t                  get_global_hid(){return global_hid++;}

  dvec&                   get_psi(){return psi;}
  dvec&                   get_EIR(){return EIR;}

  void                    zero_CC(){CC = 0.;}
  void                    inc_CC(double CC_delta){CC += CC_delta;}
  double                  get_CC(){return CC;}

  void                    zero_WW(){WW = 0.;}
  void                    inc_WW(double WW_delta){WW += WW_delta;}
  double                  get_WW(){return WW;}

  void                    zero_ZZ(){ZZ = 0.;}
  void                    inc_ZZ(double ZZ_delta){ZZ += ZZ_delta;}
  double                  get_ZZ(){return ZZ;}

  Rcpp::List              get_state_age(){return state_age;}

  Rcpp::IntegerMatrix&    get_state_age_tnow(){return state_age_tnow;}

  Rcpp::IntegerMatrix&    get_cinc_age(){return cinc_age;}
  void                    push_cinc_age(const size_t j);

  Rcpp::IntegerMatrix&    get_mosquito(){return mosquito;}
  void                    push_mosy(const int S, const int E, const int I);

  dvec&                   get_lambda_v(){return lambda_v;}
  void                    push_lambda_v(const double lambda_v_t);


private:

  /* constructor & destructor */
  globals();
  ~globals();

  /* delete all copy & move semantics */
  globals(const globals&) = delete;
  globals& operator=(const globals&) = delete;
  globals(globals&&) = delete;
  globals& operator=(globals&&) = delete;

  size_t                  NHOUSE;

  // hold the things that calculate means/variances
  stat_map                running_stats;

  // parameters
  pmap                    parameters;

  /* state (required for calculation of hazards) */
  size_t                  tnow;
  size_t                  tmax;
  size_t                  global_hid;

  // MOSQUITO -> HUMAN
  dvec                    psi;
  dvec                    EIR;

  // HUMAN -> MOSQUITO
  double                  CC; // P(bite will cause infection in mosquito) --- expectation of this prob over all landscape/individual heterogeneities
  double                  WW; // avg probability to bite and survive
  double                  ZZ; // avg probability to bite

  /* output */
  // state
  Rcpp::List              state_age; // (states X age) by time
  Rcpp::IntegerMatrix     state_age_tnow; // states X age (for this time-step, push it into the list 'state_age' before we leave this iteration)
  Rcpp::IntegerMatrix     cinc_age; // clinical incidence X age by time
  Rcpp::IntegerMatrix     mosquito; // SEI

  // rates
  dvec                    lambda_v;
  Rcpp::NumericMatrix     eir;
  Rcpp::NumericMatrix     b;

};

#endif
