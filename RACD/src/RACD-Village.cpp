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
 #  Village Class Implementation
*/

#include "RACD-Village.hpp"
#include "RACD-House.hpp"
#include "RACD-Human.hpp"
#include "RACD-Parameters.hpp"
#include "RACD-PRNG.hpp"


/* iterative mean: Knuth, The Art of Computer Programming Vol 2, section 4.2.2 */
inline double mean(const std::vector<double> array){
  double avg = 0;
  int t = 1;
  for (double x : array) {
    avg += (x - avg) / t;
    ++t;
  }
  return avg;
}

/* constructor */
village::village(const Rcpp::List& human_par, const Rcpp::List& house_par){

  #ifdef DEBUG_HPP
  std::cout << "village being born at " << this << std::endl;
  #endif

  /* houses */
  for(size_t i=0; i<house_par.size(); i++){
    Rcpp::List hh = house_par[i];
    houses.push_back(std::make_unique<house>(int(i),double(hh["psi"]),double(hh["x"]),double(hh["y"])));
  }

  /* humans */
  for(size_t i=0; i<human_par.size(); i++){

    /* human i's parameters */
    Rcpp::List hh = human_par[i];
    int humanID = int(i);
    double age = double(hh["age"]);
    bool alive = bool(hh["alive"]);
    int house = int(hh["house"]) - 1; /* R is 1-indexed */
    double bitingHet = double(hh["bitingHet"]);
    double IB = double(hh["IB"]);
    double ID = double(hh["ID"]);
    double ICA = double(hh["ICA"]);
    double ICM = double(hh["ICM"]);
    double epsilon = double(hh["epsilon"]);
    double lambda = double(hh["lambda"]);
    double phi = double(hh["phi"]);
    double prDetectAMic = double(hh["prDetectAMic"]);
    double prDetectAPCR = double(hh["prDetectAPCR"]);
    double prDetectUPCR = double(hh["prDetectUPCR"]);
    std::string state(1,hh["state"]);
    int daysLatent = int(hh["daysLatent"]);

    /* add human i to their house */
    houses[house]->add_human(std::make_unique<human>(humanID,age,alive,state,daysLatent,IB,ID,ICA,ICM,bitingHet,epsilon,lambda,phi,prDetectAMic,prDetectAPCR,prDetectUPCR,houses[house].get()));
  }

};

/* destructor */
village::~village(){
  #ifdef DEBUG_HPP
  std::cout << "village being killed at " << this << std::endl;
  #endif
};

// /* Simulation Methods */
//
// /* daily simulation */
// void                                      one_day();
//
/* demographics */
void village::births(){

  /* current population size */
  int N = 0;
  for(auto &h : houses){
    N += h->get_humans().size();
  }

  /* number of births */
  double mu = RACD_Parameters::instance()->get_mu();
  int numNewBirths = RACD_Parameters::instance()->get_prng()->get_rbinom(N,mu);

};
