#include "globals.hpp"
#include "house.hpp"
#include "mosquito.hpp"
#include "human.hpp"

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
#include <progress.hpp>

extern Rcpp::NumericMatrix output1;

// [[Rcpp::export]]
Rcpp::NumericMatrix test(){
  output1.fill(0);
  output1.at(0,0) += R::rpois(42);
  return output1;
}


// IF WEIRD BUGS OCCUR CHANGE AROUND THE ORDER OF PLUGINS



// simulation will look something like this
//
// the daily update:
//
// call update_biting(houses) to update CC,WW,ZZ (mosquitoes need it)
// call feeding_cycle(mosy) to update EIR on houses (this uses CC,WW,ZZ to calc the FOI on mosy)
//        at the end of feeding_cycle function, humans and mosquitoes are conditionally independent of each other for this time step
// call euler_step(mosy) to run the mosy model
// call what ever update humans
// do births
// do deaths
// update household lvl interventions
// repeat until tnow > tmax


// [[Rcpp::export]]
void tiny_racd(
  const Rcpp::NumericVector& theta,
  const size_t tmax
){

  /* clear global variables */
  reset_globals(tmax);

  /* put parameters in hash table */
  Rcpp::CharacterVector theta_names = theta.names();
  for(size_t i=0; i<theta.size(); i++){
    parameters.emplace(theta_names.at(i),theta.at(i));
  }

};
