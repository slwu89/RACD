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
inline double iter_mean(const std::vector<double>& array){
  double avg = 0;
  int t = 1;
  for (double x : array) {
    avg += (x - avg) / t;
    ++t;
  }
  return avg;
}

/* comparator to check if the human is dead */
bool died(const human& h){
  return !h.get_alive();
};


/* constructor */
village::village(const Rcpp::List& human_par, const Rcpp::List& house_par){

  #ifdef DEBUG_HPP
  std::cout << "village being born at " << this << std::endl;
  #endif

  max_humanID = 0;

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

    max_humanID++;
  }

};

/* destructor */
village::~village(){
  #ifdef DEBUG_HPP
  std::cout << "village being killed at " << this << std::endl;
  #endif
};

/* Simulation Methods */

/* one simulation run */
void village::simulation(const int& tMax){

  Progress::Progress pb(tMax,true);
std::cout << "entering simulation..." << std::endl;
  for(int i=0; i<tMax; i++){
    std::cout << "i: " << i << std::endl;
    if(Progress::check_abort()){
      Rcpp::stop("user abort detected; exiting RACD");
    }
    one_day();
    pb.increment();
  }

};

/* daily simulation */
void village::one_day(){

  /* run daily simulation for all humans */
  for(auto &hh : houses){
    for(auto &h : hh->get_humans()){
      h->one_day();
    }
  }

  /* demographics */
  births();
  deaths();

};

/* demographics */
void village::births(){

  /* current population size */
  int N = 0;
  for(auto &hh : houses){
    N += hh->get_humans().size();
  }
  std::cout << "got current pop: " << N << std::endl;
  /* number of births */
  double mu = RACD_Parameters::instance()->get_mu();
  int numNewBirths = RACD_Parameters::instance()->get_prng()->get_rbinom(N,mu);
std::cout << "got num births: " << numNewBirths << std::endl;
  /* ICM for newborns is function of ICA levels of 18-22 yr old women */
  std::vector<double> ICA18_22;
  for(auto &hh : houses){
    for(auto &h : hh->get_humans()){
      double age = h->get_age();
      if(age >= 18 && age < 22){
        ICA18_22.push_back(h->get_ICA());
      }
    }
  }

  double meanICA18_22 = iter_mean(ICA18_22);
  std::cout << "got meanICA18_22: " << meanICA18_22 << std::endl;

  if(numNewBirths>0){
    /* assign newborns to smallest houses */
    for(size_t i=0; i<numNewBirths; i++){
        std::cout << "making human: " << max_humanID << std::endl;
      /* find smallest house */
      std::vector<int> hhSize;
      for(auto &hh : houses){
        hhSize.push_back(hh->get_humans().size());
      }
      auto it = std::min_element(hhSize.begin(), hhSize.end());
      size_t hh_ix = std::distance(hhSize.begin(), it);

      /* fixed parameters */
      double sigma2 = RACD_Parameters::instance()->get_sigma2();
      double PM = RACD_Parameters::instance()->get_PM();
      double epsilon0 = RACD_Parameters::instance()->get_epsilon0();
      double rho = RACD_Parameters::instance()->get_rho();
      double psi = houses[hh_ix]->get_psi();
      double b0 = RACD_Parameters::instance()->get_b0();
      double phi0 = RACD_Parameters::instance()->get_phi0();
      double phi1 = RACD_Parameters::instance()->get_phi1();
      double IC0 = RACD_Parameters::instance()->get_IC0();
      double kappaC = RACD_Parameters::instance()->get_kappaC();

      /* human i's parameters */
      int humanID = max_humanID;
      double age = 0;
      bool alive = 1;
      /* house in hh_ix */
      double bitingHet = RACD_Parameters::instance()->get_prng()->get_rlnorm(-sigma2/2,sqrt(sigma2));
      double IB = 0;
      double ID = 0;
      double ICA = 0;
      double ICM = PM * meanICA18_22;
      double epsilon = epsilon0 * bitingHet * (1-rho) * psi;
      double lambda = epsilon * b0;
      double phi = phi0 * (phi1 + ((1 - phi1)/(1 + pow((PM*meanICA18_22/IC0),kappaC))));
      double prDetectAMic = 1;
      double prDetectAPCR = 1;
      double prDetectUPCR = 1;
      std::string state = "S";
      int daysLatent = 0;

      /* add human i to their house */
      houses[hh_ix]->add_human(std::make_unique<human>(humanID,age,alive,state,daysLatent,IB,ID,ICA,ICM,bitingHet,epsilon,lambda,phi,prDetectAMic,prDetectAPCR,prDetectUPCR,houses[hh_ix].get()));

      /* increment maximum humanID */
      max_humanID++;
    } /* end for */
  } /* end if */

};

/* clear out corpses at end of each day */
void village::deaths(){

  for(auto &hh : houses){

    auto it = std::find_if(hh->get_humans().begin(),hh->get_humans().end(), [&](human_ptr& h){ return !h->get_alive(); });
    if(it!=hh->get_humans().end()){
      hh->get_humans().erase(it);
    }
    it++;

  }

};
