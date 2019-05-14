#include "globals.hpp"
#include "house.hpp"
#include "mosquito.hpp"
#include "human.hpp"

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppProgress)]]

#include <memory>
#include <string>
#include <vector>

#include <Rcpp.h>
#include <progress.hpp>

// the houses (houses is static btwn calls)
using house_ptr = std::unique_ptr<house>;
using house_vector = std::vector<house_ptr>;
static house_vector houses;


// the mosquitos
using mosquito_ptr = std::unique_ptr<mosquitos>;


// IF WEIRD BUGS OCCUR CHANGE AROUND THE ORDER OF PLUGINS

// /* shamelessly "referenced" from Knuth, The Art of Computer Programming Vol 2, section 4.2.2 */
// // the mean immunity among 18-22 year olds
// inline double mean_ICA18_22(){
//   double avg = 0.;
//   int t = 1;
//   for (auto& hh : houses) {
//     for(auto& h : hh->humans){
//       if((h->age >= 18.) && (h->age < 22.)){
//         avg += (h->ICA - avg) / t;
//         ++t;
//       }
//     }
//   }
//   // weird hack if there aren't any 18-22 yr olds
//   // just take avg ICA of entire pop
//   if(avg < 0.001){
//     avg = 0.;
//     t = 1;
//     for (auto& hh : houses) {
//       for(auto& h : hh->humans){
//         avg += (h->ICA - avg) / t;
//         ++t;
//       }
//     }
//   }
//   return avg;
// }


// population level tracking

// // track states
// inline void track_state(){
//
//   for(auto& hh : houses){
//     for(auto& h : hh->humans){
//       if(h->state.compare("S") == 0){
//         state_S.at(tnow) += 1;
//       } else if(h->state.compare("E") == 0){
//         state_E.at(tnow) += 1;
//       } else if(h->state.compare("T") == 0){
//         state_T.at(tnow) += 1;
//       } else if(h->state.compare("D") == 0){
//         state_D.at(tnow) += 1;
//       } else if(h->state.compare("A") == 0){
//         state_A.at(tnow) += 1;
//       } else if(h->state.compare("U") == 0){
//         state_U.at(tnow) += 1;
//       } else if(h->state.compare("P") == 0){
//         state_P.at(tnow) += 1;
//       } else {
//         Rcpp::stop("incorrect state detected");
//       }
//     }
//   }
//
// }
//
// // track age group
// inline void track_age(){
//
//   for(auto& hh : houses){
//     for(auto& h : hh->humans){
//
//       if((h->age >= 2.) && (h->age < 10.)){
//         num_2_10.at(tnow) += 1;
//       }
//       if(h->age < 5.) {
//         num_0_5.at(tnow) += 1;
//       } else if((h->age >= 5.) && (h->age < 10.)){
//         num_5_10.at(tnow) += 1;
//       } else if((h->age >= 10.) && (h->age < 15.)){
//         num_10_15.at(tnow) += 1;
//       } else if(h->age >= 15.){
//         num_15Plus.at(tnow) += 1;
//       }
//
//       num_All.at(tnow) += 1;
//
//     }
//   }
//
// };


// // daily births
// inline void one_day_births(){
//
//   size_t hpop = num_All.at(tnow);
//   double mu = parameters.at("mu");
//
//   size_t nbirth = (size_t)R::rbinom((double)hpop, mu);
//
//   if(nbirth > 0){
//
//     double ICA18_22 = mean_ICA18_22();
//
//     double sigma2 = parameters.at("sigma2");
//
//     double PM = parameters.at("PM");
//
//     double phi0 = parameters.at("phi0");
//     double phi1 = parameters.at("phi1");
//     double kappaC = parameters.at("kappaC");
//     double IC0 = parameters.at("IC0");
//
//     for(size_t i=0; i<nbirth; i++){
//
//       // put newborns in the smallest houses for ... reasons
//       size_t smallest_hh = std::distance(houses.begin(),std::min_element(houses.begin(), houses.end(), [](auto& hh1, auto& hh2) {
//         return hh1->n < hh2->n;
//       }));
//
//       // sample this person's biting heterogeneity
//       double zeta = R::rlnorm(-sigma2/2., std::sqrt(sigma2));
//
//       // P(clinical disease | infection)
//       double phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(PM*ICA18_22/IC0,kappaC))));
//
//       // put the human in their new home
//       houses.at(smallest_hh)->humans.emplace_back(std::make_unique<human>(
//         0.,
//         houses.at(smallest_hh).get(),
//         zeta,
//         0.,
//         0.,
//         0.,
//         (PM * ICA18_22),
//         phi,
//         1.,
//         1.,
//         1.,
//         "S"
//       ));
//
//     }
//
//   }
//
// }


// simulation will look something like this
//
// the daily update:
//
// call update_biting(houses) to update CC,WW,ZZ (mosquitoes need it)
// call feeding_cycle(mosy) to update EIR on houses (this uses CC,WW,ZZ to calc the FOI on mosy)
//        at the end of feeding_cycle function, humans and mosquitoes are conditionally independent of each other for this time step
// call euler_step(mosy) to run the mosy model
// call what ever update humans
// do births
// do deaths
// update household lvl interventions
// repeat until tnow > tmax





// [[Rcpp::export]]
void tiny_racd(
  const Rcpp::List& humans_param,
  const Rcpp::List& house_param,
  const Rcpp::List& mosy_param,
  const Rcpp::NumericVector& theta,
  const size_t tmax
){

  Rcpp::Rcout << " --- initializing global variables and parameters --- " << std::endl;

  /* clear global variables */
  reset_globals(tmax);

  houses.clear();

  /* put parameters in hash table */
  Rcpp::CharacterVector theta_names = theta.names();
  for(size_t i=0; i<theta.size(); i++){
    parameters.emplace(theta_names.at(i),theta.at(i));
  }

  Rcpp::Rcout << " --- done initializing global variables and parameters --- " << std::endl;

  Rcpp::Rcout << " --- initializing simulation objects in memory --- " << std::endl;

  /* make the houses */
  size_t nhouse = house_param.size();
  psi.reserve(nhouse);
  EIR.reserve(nhouse);
  for(size_t i=0; i<nhouse; i++){

    /* biting weight */
    double psi_i = Rcpp::as<double>(Rcpp::as<Rcpp::List>(house_param[i])["psi"]);
    psi.emplace_back(psi_i);

    /* make the house */
    houses.emplace_back(
      std::make_unique<house>(i)
    );
  }

  /* make the mosquitos */
  mosquito_ptr mosy_pop = std::make_unique<mosquitos>(
    Rcpp::as<int>(mosy_param["EL_eq"]),
    Rcpp::as<int>(mosy_param["LL_eq"]),
    Rcpp::as<int>(mosy_param["PL_eq"]),
    Rcpp::as<int>(mosy_param["SV_eq"]),
    Rcpp::as<int>(mosy_param["EV_eq"]),
    Rcpp::as<int>(mosy_param["IV_eq"]),
    Rcpp::as<int>(mosy_param["K_eq"])
  );

  /* make the humans */
  size_t nhum = humans_param.size();
  for(size_t i=0; i<nhum; i++){
    Rcpp::List hum_i = Rcpp::as<Rcpp::List>(humans_param[i]);
    size_t house_id = Rcpp::as<size_t>(hum_i["house"]) - 1;
    houses[house_id]->humans.emplace_back(
      std::make_unique<human>(
        Rcpp::as<double>(hum_i["age"]),
        houses[house_id].get(),
        Rcpp::as<double>(hum_i["zeta"]),
        Rcpp::as<double>(hum_i["IB"]),
        Rcpp::as<double>(hum_i["ID"]),
        Rcpp::as<double>(hum_i["ICA"]),
        Rcpp::as<double>(hum_i["ICM"]),
        Rcpp::as<double>(hum_i["phi"]),
        Rcpp::as<double>(hum_i["prDetectAMic"]),
        Rcpp::as<double>(hum_i["prDetectAPCR"]),
        Rcpp::as<double>(hum_i["prDetectUPCR"]),
        Rcpp::as<std::string>(hum_i["state"])
      )
    );
  }

  Rcpp::Rcout << " --- done initializing simulation objects in memory --- " << std::endl;


  Rcpp::Rcout << " --- begin simulation --- " << std::endl;

  // main simulation loop
  Progress pb(tmax,true);
  while(tnow < tmax){

    // check for worried users
    if(tnow % 5 == 0){
      if(Progress::check_abort()){
        Rcpp::stop("user abort detected");
      }
    }

    // track human output
    track_state(houses);
    track_age(houses);

    // track mosquito output
    track_mosquito(mosy_pop);

    // update house -> mosy state variables
    update_biting(houses);

    // run mosquito biting (mosy -> house transmission)
    feeding_cycle(mosy_pop);

    // AFTER THIS POINT HUMANS/MOSY ARE CONDITIONALLY INDEPENDENT OF EACH OTHER

    // mosquito sim
    euler_step(mosy_pop);


    // bookkeeping before we move on
    pb.increment();
    tnow++;
  }


  Rcpp::Rcout << std::endl << " --- end simulation --- " << std::endl;


};
