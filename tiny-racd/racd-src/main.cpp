#include "globals.hpp"

extern Rcpp::NumericMatrix output1;

// [[Rcpp::export]]
Rcpp::NumericMatrix test(){
  output1.fill(0);
  output1.at(0,0) += R::rpois(42);
  return output1;
}




// simulation will look something like this
//
// the daily update:
//
// call update_biting(houses) to update CC,WW,ZZ (mosquitoes need it)
// call feeding_cycle(mosy) to update EIR on houses (this uses CC,WW,ZZ to calc the FOI on mosy)
// call euler_step(mosy) to run the mosy model
// call what ever update humans
// do births
// do deaths
// update household lvl interventions
// repeat until tnow > tmax
