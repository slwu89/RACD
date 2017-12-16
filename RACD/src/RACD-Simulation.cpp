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

//' RACD Simulation Model
//'
//' Run the simulation
//'
//' @param tMax length of simulation in days
//' @param theta named list of parameters (see \code{\link[RACDaux]{RACD_Parameters}} for details)
//' @param human list of human parameters (see \code{\link[RACDaux]{RACD_Setup}} for details)
//' @param house list of house parameters (see \code{\link[RACDaux]{RACD_Setup}} for details)
//' @param seed seed for prng class
//'
//' @examples
//' theta = RACDaux::RACD_Parameters()
//' init = RACDaux::RACD_Setup(theta)
//' \dontrun{RACD_Simulation(365,theta,init$humans,init$houses,123)}
//'
//' @export
// [[Rcpp::export]]
void RACD_Simulation(const int& tMax, const Rcpp::NumericVector &theta, const Rcpp::List& human, const Rcpp::List& house, const uint_least32_t &seed){

  /* initialize parameters and RNG seed */
  RACD_Parameters::instance();
  RACD_Parameters::instance()->set_values(theta["epsilon0"],theta["fT"],theta["dE"],theta["dT"],theta["dD"],theta["dA"],theta["dU"],theta["dP"],theta["cD"],theta["cT"],theta["cU"],theta["gammaI"],theta["rho"],theta["a0"],theta["sigma2"],theta["d1"],theta["dID"],theta["ID0"],theta["kappaD"],theta["uD"],theta["aD"],theta["fD0"],theta["gammaD"],theta["alphaA"],theta["alphaU"],theta["b0"],theta["b1"],theta["dB"],theta["IB0"],theta["kappaB"],theta["uB"],theta["phi0"],theta["phi1"],theta["dC"],theta["IC0"],theta["kappaC"],theta["uC"],theta["PM"],theta["dM"],theta["rW"],theta["rP"],theta["meanAge"],theta["N"],theta["meanNumPeoplePerHouse"],theta["numHousesPerBreedingSite"]);
  RACD_Parameters::instance()->set_prng(seed);

  /* create village */
  village simVillage(human,house);
  simVillage.simulation(tMax);

  /* kill global singleton */
  RACD_Parameters::instance()->suicide();
};
