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
 #  R Interface to RACD Simulation
*/

/* C++ includes */
#include <memory>

/* Rcpp includes */
#include <Rcpp.h>

#include <progress.hpp>
#include <progress_bar.hpp>

/* RACD includes */
#include "RACD-Human.hpp"
#include "RACD-House.hpp"
#include "RACD-Village.hpp"
#include "RACD-Mosquito.hpp"

//' RACD Simulation Model
//'
//' Run the simulation
//'
//' @param tMax length of simulation in days
//' @param theta named list of parameters (see \code{\link{RACD_Parameters}} for details)
//' @param human list of human parameters (see \code{\link{RACD_Setup}} for details)
//' @param house list of house parameters (see \code{\link{RACD_Setup}} for details)
//' @param seed seed for prng class
//' @param outfile path to .csv file for logging state transition events
//'
//' @examples
//' \dontrun{
//' library(RACD)
//' library(tidyverse)
//' library(spatstat)
//' xy_h <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
//' xy_b <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
//' theta <- RACD_Parameters()
//' init <- RACD_Setup(as.matrix.ppx(xy_h),as.matrix.ppx(xy_b),theta)
//' outfile = "/Users/slwu89/Desktop/log_trans.csv"
//' RACD_Simulation(365,theta,init$humans,init$houses,123,outfile)
//' state = RACD_StateVector(outfile)
//' state %>% as.tibble %>% gather(state,value,-time) %>% ggplot(aes(x=time,y=value,color=state)) + geom_line() + theme_bw()
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List RACD_Simulation(const int tMax, const Rcpp::List &theta, const Rcpp::List& human, const Rcpp::IntegerVector& mosy, const Rcpp::List& house){

  /* output */
  Rcpp::IntegerMatrix state_out(tMax,6);
  Rcpp::colnames(state_out) = Rcpp::CharacterVector::create("S","T","D","A","U","P");

  Rcpp::IntegerMatrix mosy_out(tMax,6);
  Rcpp::colnames(mosy_out) = Rcpp::CharacterVector::create("EL","LL","PL","SV","EV","IV");

  /* construct village */
  std::unique_ptr<village> village_ptr(std::make_unique<village>(theta,mosy));

  /* initialize objects */
  village_ptr->initialize(human,house);

  /* run simulation */
  Progress pb(tMax,true);
  while(village_ptr->tNow < tMax){

    if(village_ptr->tNow % 2 == 0){
      Rcpp::checkUserInterrupt();
    }

    village_ptr->one_day();

    /* track output */
    village_ptr->track_human_state(state_out);
    village_ptr->track_mosquito_state(mosy_out);

    /* increment time */
    pb.increment();
    village_ptr->tNow++;
  }

  return Rcpp::List::create(
    Rcpp::Named("H_state") = state_out,
    Rcpp::Named("M_state") = mosy_out
  );
};
