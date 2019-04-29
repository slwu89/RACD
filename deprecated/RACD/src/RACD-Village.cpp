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


/* ######################################################################
 # RACD includes & auxiliary functions
###################################################################### */

#include "RACD-Village.hpp"
#include "RACD-House.hpp"
#include "RACD-Human.hpp"
#include "RACD-Parameters.hpp"
#include "RACD-PRNG.hpp"
// #include "RACD-Logger.hpp"


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

/* lookup table for mapping human states to columns of output matrix */
static const std::unordered_map<std::string,size_t> state2col = {
  {"S",0},
  {"E",1},
  {"T",2},
  {"D",3},
  {"A",4},
  {"U",5},
  {"P",6},
};


/* ######################################################################
 # construtor & destructor
###################################################################### */

/* constructor */
village::village(const uint_least32_t seed, const Rcpp::NumericVector& theta) :

  /* initialize utility classes */
  prng_ptr(std::make_unique<prng>(seed)),
  // logger_ptr(std::make_unique<logger>()),
  param_ptr(std::make_unique<parameters>(
    Rcpp::as<double>(theta["epsilon0"]),
    Rcpp::as<double>(theta["fT"]),
    Rcpp::as<int>(theta["dE"]),
    Rcpp::as<int>(theta["dT"]),
    Rcpp::as<int>(theta["dD"]),
    Rcpp::as<int>(theta["dA"]),
    Rcpp::as<int>(theta["dU"]),
    Rcpp::as<int>(theta["dP"]),
    Rcpp::as<double>(theta["cD"]),
    Rcpp::as<double>(theta["cT"]),
    Rcpp::as<double>(theta["cU"]),
    Rcpp::as<double>(theta["gammaI"]),
    Rcpp::as<double>(theta["rho"]),
    Rcpp::as<double>(theta["a0"]),
    Rcpp::as<double>(theta["sigma2"]),
    Rcpp::as<double>(theta["d1"]),
    Rcpp::as<double>(theta["dID"]),
    Rcpp::as<double>(theta["ID0"]),
    Rcpp::as<double>(theta["kappaD"]),
    Rcpp::as<double>(theta["uD"]),
    Rcpp::as<double>(theta["aD"]),
    Rcpp::as<double>(theta["fD0"]),
    Rcpp::as<double>(theta["gammaD"]),
    Rcpp::as<double>(theta["alphaA"]),
    Rcpp::as<double>(theta["alphaU"]),
    Rcpp::as<double>(theta["b0"]),
    Rcpp::as<double>(theta["b1"]),
    Rcpp::as<double>(theta["dB"]),
    Rcpp::as<double>(theta["IB0"]),
    Rcpp::as<double>(theta["kappaB"]),
    Rcpp::as<double>(theta["uB"]),
    Rcpp::as<double>(theta["phi0"]),
    Rcpp::as<double>(theta["phi1"]),
    Rcpp::as<double>(theta["dC"]),
    Rcpp::as<double>(theta["IC0"]),
    Rcpp::as<double>(theta["kappaC"]),
    Rcpp::as<double>(theta["uC"]),
    Rcpp::as<double>(theta["PM"]),
    Rcpp::as<double>(theta["dM"]),
    Rcpp::as<double>(theta["rW"]),
    Rcpp::as<double>(theta["rP"]),
    Rcpp::as<double>(theta["meanAge"]),
    Rcpp::as<int>(theta["N"]))
  ),

  /* initialize data members */
  tNow(0),
  run_id(1),
  max_humanID(0)
{
  #ifdef DEBUG_RACD
  std::cout << "village ctor being called at " << this << std::endl;
  #endif
};


/* destructor */
village::~village(){
  #ifdef DEBUG_RACD
  std::cout << "village dtor being called at " << this << std::endl;
  #endif
};


/* ######################################################################
 # initialize objects
###################################################################### */

void village::initialize(const Rcpp::List &humansR, const Rcpp::List &housesR){

  Rcpp::Rcout << "initializing simulation objects ... ";

  /* houses */
  for(size_t i=0; i<housesR.size(); i++){
    houses.emplace_back(std::make_unique<house>(
      int(i),
      Rcpp::as<double>(Rcpp::as<Rcpp::List>(housesR[i])["psi"]),
      Rcpp::as<double>(Rcpp::as<Rcpp::List>(housesR[i])["x"]),
      Rcpp::as<double>(Rcpp::as<Rcpp::List>(housesR[i])["y"]),
      this
    ));
  }

  /* humans */
  for(size_t i=0; i<humansR.size(); i++){

    /* parameters needed to log and place the human */
    int house = Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["house"]) - 1;
    int id = max_humanID;
    double age = Rcpp::as<double>(Rcpp::as<Rcpp::List>(humansR[i])["age"]);
    std::string state = Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(humansR[i])["state"]);

    /* log human's initial state */
    // logger_ptr->get_log() << std::to_string(id) << "," << state << ",0," << std::to_string(age) << "\n";

    /* put them in their house */
    houses[house]->add_human(std::make_unique<human>(
      id,
      age,
      Rcpp::as<bool>(Rcpp::as<Rcpp::List>(humansR[i])["alive"]),
      state,
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["daysLatent"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["IB"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["ID"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["ICA"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["ICM"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["bitingHet"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["epsilon"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["lambda"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["phi"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["prDetectAMic"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["prDetectAPCR"]),
      Rcpp::as<int>(Rcpp::as<Rcpp::List>(humansR[i])["prDetectUPCR"]),
      houses[house].get())
    );

    /* increment total number of human id */
    max_humanID++;
  }

  Rcpp::Rcout << "done!" << std::endl;

};


/* ######################################################################
 # Simulation methods
###################################################################### */

/* one simulation run */
void village::simulation(const int tMax, Rcpp::IntegerMatrix& state_out){

  Progress pb(tMax,true);

  /* run simulation */
  while(tNow < tMax){

    // output
    for(auto &hh : houses){
      for(auto &h : hh->get_humans()){
        size_t col = state2col.at(h->get_state());
        state_out.at(tNow, col) += 1;
      }
    }

    if(tNow % 20 == 0){
      Rcpp::checkUserInterrupt();
    }

    one_day();
    pb.increment();

    tNow++;
  }

};


/* daily simulation */
void village::one_day(){

  /* run daily simulation for all humans */
  for(auto &hh : houses){
    for(auto &h : hh->get_humans()){
      h->one_day(tNow);
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

  /* number of births */
  double mu = param_ptr->get_mu();
  int numNewBirths = prng_ptr->get_rbinom(N,mu);

  if(numNewBirths>0){

    /* ICM for newborns is function of ICA levels of 18-22 yr old women */
    double meanICA18_22 = 0.0;
    int t = 1;
    // std::vector<double> ICA18_22;
    for(auto &hh : houses){
      for(auto &h : hh->get_humans()){
        double age = h->get_age();
        if(age >= 18 && age < 22){
          // ICA18_22.push_back(h->get_ICA());
          /* iterative mean: Knuth, The Art of Computer Programming Vol 2, section 4.2.2 */
          meanICA18_22 += (h->get_ICA() - meanICA18_22) / t;
          ++t;
        }
      }
    }

    // double meanICA18_22 = iter_mean(ICA18_22);

    /* assign newborns to smallest houses */
    for(size_t i=0; i<numNewBirths; i++){

      /* find smallest house */
      std::vector<int> hhSize;
      for(auto &hh : houses){
        hhSize.push_back(hh->get_humans().size());
      }
      auto it = std::min_element(hhSize.begin(), hhSize.end());
      size_t hh_ix = std::distance(hhSize.begin(), it);

      /* fixed parameters */
      double sigma2 = param_ptr->get_sigma2();
      double PM = param_ptr->get_PM();
      double epsilon0 = param_ptr->get_epsilon0();
      double rho = param_ptr->get_rho();
      double psi = houses[hh_ix]->get_psi();
      double b0 = param_ptr->get_b0();
      double phi0 = param_ptr->get_phi0();
      double phi1 = param_ptr->get_phi1();
      double IC0 = param_ptr->get_IC0();
      double kappaC = param_ptr->get_kappaC();

      /* human i's parameters */
      int humanID = max_humanID;
      double age = 0;
      bool alive = 1;
      /* house in hh_ix */
      double bitingHet = prng_ptr->get_rlnorm(-sigma2/2,sqrt(sigma2));
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

      /* logging */
      // logger_ptr->get_log() << std::to_string(humanID) << ",Birth," << std::to_string(tNow) << "," + std::to_string(age) << "\n";

      /* add human i to their house */
      houses[hh_ix]->add_human(std::make_unique<human>(humanID,age,alive,state,daysLatent,IB,ID,ICA,ICM,bitingHet,epsilon,lambda,phi,prDetectAMic,prDetectAPCR,prDetectUPCR,houses[hh_ix].get()));

      /* increment maximum humanID */
      max_humanID++;
    } /* end for */
  } /* end if */

};


/* clear out dead humans at end of each day */
void village::deaths(){

  for(auto &hh : houses){

    auto it = std::find_if(hh->get_humans().begin(),hh->get_humans().end(), [&](human_ptr& h){ return !h->get_alive(); });
    if(it!=hh->get_humans().end()){
      hh->get_humans().erase(it);
    }
    it++;

  }

};
