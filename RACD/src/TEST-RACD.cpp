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
 #  Anomalous Materials Lab
*/

#include <Rcpp.h>
#include "RACD-Parameters.hpp"
#include "RACD-PRNG.hpp"
#include "TEST-HUMAN.hpp"

#include <iostream>


//' Unit Test: RACD_Parameters
//'
//' @param theta output of \code{\link[RACDaux]{RACD_Parameters}}
//'
//' @export
// [[Rcpp::export]]
void test_parameters(const Rcpp::NumericVector &theta){
  std::cout << "begin testing singleton parameters" << std::endl;
  RACD_Parameters::instance();
  RACD_Parameters::instance()->set_values(
    theta["epsilon0"],
    theta["fT"],
     theta["dE"],
     theta["dT"],
     theta["dD"],
     theta["dA"],
     theta["dU"],
     theta["dP"],
     theta["cD"],
     theta["cT"],
     theta["cU"],
     theta["gammaI"],
     theta["rho"],
     theta["a0"],
     theta["sigma2"],
     theta["d1"],
     theta["dID"],
     theta["ID0"],
     theta["kappaD"],
     theta["uD"],
     theta["aD"],
     theta["fD0"],
     theta["gammaD"],
     theta["alphaA"],
     theta["alphaU"],
     theta["b0"],
     theta["b1"],
     theta["dB"],
     theta["IB0"],
     theta["kappaB"],
     theta["uB"],
     theta["phi0"],
     theta["phi1"],
     theta["dC"],
     theta["IC0"],
     theta["kappaC"],
     theta["uC"],
     theta["PM"],
     theta["dM"],
     theta["rW"],
     theta["rP"],
     theta["meanAge"],
     theta["N"],
     theta["meanNumPeoplePerHouse"],
     theta["numHousesPerBreedingSite"]
  );
  std::cout << "test getting kappaD: " << RACD_Parameters::instance()->get_kappaD() << std::endl;
  RACD_Parameters::instance()->suicide();
  std::cout << "exit testing singleton parameters" << std::endl;
}


//' Unit Test: test_human
//'
//' @param theta output of \code{\link[RACDaux]{RACD_Parameters}}
//'
//' @export
// [[Rcpp::export]]
void test_human_parameters(const Rcpp::NumericVector &theta){
  std::cout << "begin testing getting singleton parameters from a test_human" << std::endl;
  RACD_Parameters::instance();
  RACD_Parameters::instance()->set_values(
    theta["epsilon0"],
    theta["fT"],
     theta["dE"],
     theta["dT"],
     theta["dD"],
     theta["dA"],
     theta["dU"],
     theta["dP"],
     theta["cD"],
     theta["cT"],
     theta["cU"],
     theta["gammaI"],
     theta["rho"],
     theta["a0"],
     theta["sigma2"],
     theta["d1"],
     theta["dID"],
     theta["ID0"],
     theta["kappaD"],
     theta["uD"],
     theta["aD"],
     theta["fD0"],
     theta["gammaD"],
     theta["alphaA"],
     theta["alphaU"],
     theta["b0"],
     theta["b1"],
     theta["dB"],
     theta["IB0"],
     theta["kappaB"],
     theta["uB"],
     theta["phi0"],
     theta["phi1"],
     theta["dC"],
     theta["IC0"],
     theta["kappaC"],
     theta["uC"],
     theta["PM"],
     theta["dM"],
     theta["rW"],
     theta["rP"],
     theta["meanAge"],
     theta["N"],
     theta["meanNumPeoplePerHouse"],
     theta["numHousesPerBreedingSite"]
  );

  std::unique_ptr<test_human> h = std::make_unique<test_human>("sam");
  h->get_kappaD();

  RACD_Parameters::instance()->suicide();
  std::cout << "exit testing getting singleton parameters from a test_human" << std::endl;
}


//' Unit Test: prng
//'
//' test the composition prng class within the parameters singleton
//'
//' @export
// [[Rcpp::export]]
void test_prng(const uint_least32_t &seed){
  std::cout << "begin testing prng composition" << std::endl;
  RACD_Parameters::instance();
  RACD_Parameters::instance()->set_prng(seed);
  std::cout << "get some runifs ";
  for(size_t i=0; i<10; i++){
    std::cout << RACD_Parameters::instance()->get_prng()->get_runif();
    std::cout <<  " <--> ";
  }
  std::cout << std::endl;
  std::cout << "exit testing prng composition" << std::endl;
};