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
#include "RACD-Human.hpp"
#include "RACD-House.hpp"
#include "RACD-Village.hpp"

#include "TEST-HUMAN.hpp" /* test human object */

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


//' Unit Test: house
//'
//' test putting a human in a house and calling some basic actions
//'
//' @examples
//' test_house(theta = RACDaux::RACD_Parameters(),seed = 1)
//' @export
// [[Rcpp::export]]
void test_house(const Rcpp::NumericVector &theta, const uint_least32_t &seed){

  std::cout << "begin testing house & human" << std::endl;
  RACD_Parameters::instance();
  RACD_Parameters::instance()->set_values(theta["epsilon0"],theta["fT"],theta["dE"],theta["dT"],theta["dD"],theta["dA"],theta["dU"],theta["dP"],theta["cD"],theta["cT"],theta["cU"],theta["gammaI"],theta["rho"],theta["a0"],theta["sigma2"],
    theta["d1"],theta["dID"],theta["ID0"],theta["kappaD"],theta["uD"],theta["aD"],theta["fD0"],theta["gammaD"],theta["alphaA"],theta["alphaU"],theta["b0"],theta["b1"],theta["dB"],theta["IB0"],theta["kappaB"],theta["uB"],theta["phi0"],theta["phi1"],theta["dC"],theta["IC0"],theta["kappaC"],theta["uC"],theta["PM"],theta["dM"],theta["rW"],theta["rP"],theta["meanAge"],theta["N"],theta["meanNumPeoplePerHouse"],theta["numHousesPerBreedingSite"]
  );
  RACD_Parameters::instance()->set_prng(seed);



  house aHouse(1,1,1,1);

  human_ptr bob = std::make_unique<human>(1,1,1,"S",1,1,1,1,1,1,1,1,1,1,1,1,&aHouse);
  human_ptr alice = std::make_unique<human>(2,1,1,"E",15,1,1,1,1,1,1,1,1,1,1,1,&aHouse);

  aHouse.add_human(std::move(bob));
  aHouse.add_human(std::move(alice));
  aHouse.add_human(std::make_unique<human>(3,1,1,"T",1,1,1,1,1,1,1,1,1,1,1,1,&aHouse));

  std::cout << "there are " << aHouse.get_humans().size() << " humans in the house" << std::endl;

  for(auto &h : aHouse.get_humans()){
    std::cout << "iterating..." << std::endl;
    h->one_day(100);
  }
  aHouse.get_humans().front()->suicide();

  std::cout << "there are " << aHouse.get_humans().size() << " humans in the house" << std::endl;



  std::cout << "exit testing house & human" << std::endl;
}

//' Unit Test: nested lists
//'
//' test iteration through nested list
//'
//' @examples
//' theta = RACDaux::RACD_Parameters()
//' init = RACDaux::RACD_Setup(theta)
//' test_list(init$humans)
//' @export
// [[Rcpp::export]]
void test_list(const Rcpp::List& list){
  for(size_t i=0; i<list.size(); i++){
    Rcpp::List h = list[i];
    std::cout << "age: " << double(h["age"]) << std::endl;
    std::cout << "alive: " << bool(h["alive"]) << std::endl;
    std::cout << "state: " << char(h["state"]) << std::endl;
    std::string s(1,h["state"]);
    std::cout << "state as string: " << s << std::endl;
  }
}


//' Unit Test: make a village
//'
//' create a village
//'
//' @examples
//' theta = RACDaux::RACD_Parameters()
//' init = RACDaux::RACD_Setup(theta)
//' human_par = init$humans
//' house_par = init$houses
//' @export
// [[Rcpp::export]]
void test_village(const Rcpp::List& human_par, const Rcpp::List& house_par, const Rcpp::NumericVector &theta, const uint_least32_t &seed){
  RACD_Parameters::instance();
  RACD_Parameters::instance()->set_values(theta["epsilon0"],theta["fT"],theta["dE"],theta["dT"],theta["dD"],theta["dA"],theta["dU"],theta["dP"],theta["cD"],theta["cT"],theta["cU"],theta["gammaI"],theta["rho"],theta["a0"],theta["sigma2"],theta["d1"],theta["dID"],theta["ID0"],theta["kappaD"],theta["uD"],theta["aD"],theta["fD0"],theta["gammaD"],theta["alphaA"],theta["alphaU"],theta["b0"],theta["b1"],theta["dB"],theta["IB0"],theta["kappaB"],theta["uB"],theta["phi0"],theta["phi1"],theta["dC"],theta["IC0"],theta["kappaC"],theta["uC"],theta["PM"],theta["dM"],theta["rW"],theta["rP"],theta["meanAge"],theta["N"],theta["meanNumPeoplePerHouse"],theta["numHousesPerBreedingSite"]);
  RACD_Parameters::instance()->set_prng(seed);

  village simVillage(human_par,house_par);
  simVillage.births();
  RACD_Parameters::instance()->suicide();
}


//' Unit Test: test rlnorm from STL <random>
//'
//' draw rlnorm
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector test_rlnorm(const double& m, const double& s, const uint_least32_t &seed){
  RACD_Parameters::instance();
  RACD_Parameters::instance()->set_prng(seed);
  Rcpp::NumericVector out;
  for(size_t i=0; i<10000; i++){
    out.push_back(RACD_Parameters::instance()->get_prng()->get_rlnorm(m,s));
  }
  return out;
}
