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

/* RACD includes */
#include "RACD-Parameters.hpp"
#include "RACD-PRNG.hpp"
#include "RACD-Human.hpp"
#include "RACD-House.hpp"
#include "RACD-Village.hpp"
#include "RACD-Logger.hpp"

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
void RACD_Simulation(const int tMax, const Rcpp::NumericVector &theta, const Rcpp::List& human, const Rcpp::List& house, const uint_least32_t seed, const std::string& outfile){

  /* construct village */
  std::unique_ptr<village> village_ptr(std::make_unique<village>(seed,theta));

  /* initialize logging */
  village_ptr->logger_ptr->open_log(outfile);
  village_ptr->logger_ptr->get_log() << "HumanID,Event,Time,Age\n";

  /* initialize objects */
  village_ptr->initialize(human,house);

  /* run simulation */
  village_ptr->simulation(tMax);

};
