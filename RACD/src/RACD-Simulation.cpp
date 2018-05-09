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

#include <Rcpp.h>
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
//' @param theta named list of parameters (see \code{\link[RACDaux]{RACD_Parameters}} for details)
//' @param human list of human parameters (see \code{\link[RACDaux]{RACD_Setup}} for details)
//' @param house list of house parameters (see \code{\link[RACDaux]{RACD_Setup}} for details)
//' @param seed seed for prng class
//' @param out_trans path to .csv file for logging state transition events
//'
//' @examples
//' \dontrun{
//' library(RACD)
//' library(RACDaux)
//' library(tidyverse)
//' library(spatstat)
//' xy_h <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
//' xy_b <- rpoispp(lambda = 100,win = owin(c(0,1),c(0,1)))
//' theta <- RACD_Parameters()
//' init <- RACD_Setup(as.matrix.ppx(xy_h),as.matrix.ppx(xy_b),theta)
//' outfile = "/Users/slwu89/Desktop/log_trans.csv"
//' RACD_Simulation(365,theta,init$humans,init$houses,123,outfile)
//' state = RACDaux::RACD_StateVector(outfile)
//' state %>% as.tibble %>% gather(state,value,-time) %>% ggplot(aes(x=time,y=value,color=state)) + geom_line() + theme_bw()
//' }
//'
//' @export
// [[Rcpp::export]]
void RACD_Simulation(const int tMax, const Rcpp::NumericVector &theta, const Rcpp::List& human, const Rcpp::List& house, const uint_least32_t seed, const std::string& out_trans){

  /* initialize parameters */
  RACD_Parameters::instance().set_values(theta["epsilon0"],theta["fT"],theta["dE"],theta["dT"],theta["dD"],theta["dA"],theta["dU"],theta["dP"],theta["cD"],theta["cT"],theta["cU"],theta["gammaI"],theta["rho"],theta["a0"],theta["sigma2"],theta["d1"],theta["dID"],theta["ID0"],theta["kappaD"],theta["uD"],theta["aD"],theta["fD0"],theta["gammaD"],theta["alphaA"],theta["alphaU"],theta["b0"],theta["b1"],theta["dB"],theta["IB0"],theta["kappaB"],theta["uB"],theta["phi0"],theta["phi1"],theta["dC"],theta["IC0"],theta["kappaC"],theta["uC"],theta["PM"],theta["dM"],theta["rW"],theta["rP"],theta["meanAge"],theta["N"],theta["meanNumPeoplePerHouse"],theta["numHousesPerBreedingSite"]);

  /* prng seed */
  prng::instance().set_seed(seed);

  /* initialize logging */
  logger::instance().open_log(out_trans);
  std::string head_trans("HumanID,Event,Time,Age");
  logger::instance().log_trans(head_trans);

  /* create village */
  village simVillage(human,house);

  /* run simulation */
  simVillage.simulation(tMax);

  /* close logging */
  logger::instance().close_log();
};
