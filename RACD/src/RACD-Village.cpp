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
#include "RACD-Mosquito.hpp"


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


/* ######################################################################
 # construtor & destructor
###################################################################### */

/* constructor */
village::village(const Rcpp::List& theta, const Rcpp::IntegerVector& mosy, const Rcpp::List& mosy_theta) :

  /* initialize data members */
  mosquito(std::make_unique<mosquito_habitat>(Rcpp::as<int>(mosy["EL"]),
                                              Rcpp::as<int>(mosy["LL"]),
                                              Rcpp::as<int>(mosy["PL"]),
                                              Rcpp::as<int>(mosy["SV"]),
                                              Rcpp::as<int>(mosy["EV"]),
                                              Rcpp::as<int>(mosy["IV"]),
                                              Rcpp::as<double>(mosy["K"]),
                                              mosy_theta)),
  tNow(0),
  max_humanID(0)
{

  /* initialize parameters */
  Rcpp::CharacterVector theta_names = theta.names();
  for(size_t i=0; i<theta.length(); i++){
    param_ptr->emplace(Rcpp::as<std::string>(theta_names[i]), Rcpp::as<double>(theta[i]));
  }

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

// /* one simulation run */
// void village::simulation(const int tMax){
//
//   Progress pb(tMax,true);
//
//   /* run simulation */
//   while(tNow < tMax){
//
//     if(tNow % 10 == 0){
//       Rcpp::checkUserInterrupt();
//     }
//
//     one_day();
//     pb.increment();
//
//     tNow++;
//   }
//
// };


/* daily simulation */
void village::one_day(){

  /* run daily simulation for all humans */
  for(auto &hh : houses){
    hh->update_intervention();
    for(auto &h : hh->get_humans()){
      h.second->one_day(tNow);
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
  double mu = param_ptr->at("mu");
  int numNewBirths = R::rbinom(N,mu);

  if(numNewBirths>0){

    /* ICM for newborns is function of ICA levels of 18-22 yr old women */
    double meanICA18_22 = 0.0;
    int t = 1;
    for(auto &hh : houses){
      for(auto &h : hh->get_humans()){
        double age = h.second->get_age();
        if(age >= 18 && age < 22){
          /* iterative mean: Knuth, The Art of Computer Programming Vol 2, section 4.2.2 */
          meanICA18_22 += (h.second->get_ICA() - meanICA18_22) / t;
          ++t;
        }
      }
    }

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
      double sigma2 = param_ptr->at("sigma2");
      double PM = param_ptr->at("PM");
      double epsilon0 = param_ptr->at("epsilon0");
      double rho = param_ptr->at("rho");
      // double psi = houses[hh_ix]->get_psi();
      double psi = 1.0;
      double b0 = param_ptr->at("b0");
      double phi0 = param_ptr->at("phi0");
      double phi1 = param_ptr->at("phi1");
      double IC0 = param_ptr->at("IC0");
      double kappaC = param_ptr->at("kappaC");

      /* human i's parameters */
      int humanID = max_humanID;
      double age = 0;
      bool alive = 1;
      /* house in hh_ix */
      double bitingHet = R::rlnorm(-sigma2/2,sqrt(sigma2));
      double IB = 0;
      double ID = 0;
      double ICA = 0;
      double ICM = PM * meanICA18_22;
      double epsilon = epsilon0 * bitingHet * (1-rho) * psi;
      double lambda = epsilon * b0;
      double phi = phi0 * (phi1 + ((1 - phi1)/(1 + std::pow((PM*meanICA18_22/IC0),kappaC))));
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
