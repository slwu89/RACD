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

#include "globals.hpp"

#include "stats.hpp"


/* constructor & destructor */
globals::globals() : NHOUSE(0), tnow(0), tmax(0), global_hid(0), CC(0.), WW(0.), ZZ(0.)
{};

globals::~globals(){};

/* utility methods */
globals& globals::instance(){
    static globals instance;
    return instance;
};


/* ################################################################################
#   utility and initialization
################################################################################ */


// set up size of output objects
void globals::set_output(const size_t tmax_){
  if(NHOUSE == 0){
    Rcpp::stop("call 'globals::set_output' after assigning 'NHOUSE'!");
  }

  tnow = 0;
  tmax = tmax_;

  // mean/var calc
  init_stats();

  // biting
  psi.assign(NHOUSE,0.);
  EIR.assign(NHOUSE,0.);

  // output: states
  state_age_tnow = Rcpp::IntegerMatrix(7,6);
  Rcpp::rownames(state_age_tnow) = Rcpp::CharacterVector::create("S","E","T","D","A","U","P");
  Rcpp::colnames(state_age_tnow) = Rcpp::CharacterVector::create("all","2-10","0-5","5-10","10-15","15+");

  state_age = Rcpp::List(tmax,state_age_tnow);

  cinc_age = Rcpp::IntegerMatrix(tmax,6);
  Rcpp::colnames(cinc_age) = Rcpp::CharacterVector::create("all","2-10","0-5","5-10","10-15","15+");

  mosquito = Rcpp::IntegerMatrix(tmax,3);
  Rcpp::colnames(mosquito) = Rcpp::CharacterVector::create("S","E","I");

  // output: rates
  lambda_v.assign(tmax,0.);

  eir = Rcpp::NumericMatrix(tmax,2);
  b = Rcpp::NumericMatrix(tmax,2);
  Rcpp::colnames(eir) = Rcpp::CharacterVector::create("mean","var");
  Rcpp::colnames(b) = Rcpp::CharacterVector::create("mean","var");
};

// call this when we're ready to return everything to R
Rcpp::List globals::get_output(){

  return Rcpp::List::create(
    Rcpp::Named("state_age") = state_age,
    Rcpp::Named("clinical_incidence") = cinc_age,
    Rcpp::Named("mosquito") = mosquito,
    Rcpp::Named("lambda_v") = lambda_v,
    Rcpp::Named("EIR") = eir,
    Rcpp::Named("b") = b
  );

};

// set up the paramters hash table
void globals::set_parameters(const Rcpp::NumericVector& theta){

  /* put parameters in hash table */
  Rcpp::CharacterVector theta_names = theta.names();

  parameters.reserve(theta_names.size());

  for(size_t i=0; i<theta.size(); i++){
    parameters.emplace(theta_names.at(i),theta.at(i));
  }
};

// running statistics
void globals::init_stats(){
  running_stats.emplace("EIR",RunningStat());
  running_stats.emplace("b",RunningStat());
};

void globals::push_stats(){

  eir.at(tnow,0) = running_stats["EIR"].Mean();
  eir.at(tnow,1) = running_stats["EIR"].Variance();

  b.at(tnow,0) = running_stats["b"].Mean();
  b.at(tnow,1) = running_stats["b"].Variance();

  running_stats["EIR"].Clear();
  running_stats["b"].Clear();
};

// called by mosquito::feeding_cycle
void globals::update_EIR(const double bites){
  for(size_t h=0; h<NHOUSE; h++){
    EIR[h] = psi[h] * bites;
  }
};

// do this at the very bottom of the loop (just before the clock ticks to tomorrow)
void globals::iterate(){
  push_stats();
  state_age.at(tnow) = state_age_tnow;
  state_age_tnow.fill(0);
  tnow++;
};


/* ################################################################################
#   logging
################################################################################ */

void globals::push_mosy(const int S, const int E, const int I){
  mosquito.at(tnow,0) = S;
  mosquito.at(tnow,1) = E;
  mosquito.at(tnow,2) = I;
};

void globals::push_lambda_v(const double lambda_v_t){
  lambda_v.at(tnow) = lambda_v_t;
};

void globals::push_cinc_age(const size_t j){
  if(j == 0 ){
    Rcpp::stop("error: 'push_cinc_age' needs j > 0 for indexing\n");
  }
  cinc_age.at(tnow,0) += 1;
  cinc_age.at(tnow,j) += 1;
};
