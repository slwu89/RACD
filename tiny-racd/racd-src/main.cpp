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
Rcpp::List tiny_racd(
  const Rcpp::List& humans_param,
  const Rcpp::List& house_param,
  const Rcpp::List& mosy_param,
  const Rcpp::NumericVector& theta,
  const size_t tmax
){

  Rcpp::Rcout << " --- initializing global variables and parameters --- " << std::endl;

  size_t nhouse = house_param.size();

  /* clear global variables */
  reset_globals(tmax,nhouse);
  houses.reserve(nhouse);

  /* put parameters in hash table */
  Rcpp::CharacterVector theta_names = theta.names();
  for(size_t i=0; i<theta.size(); i++){
    parameters.emplace(theta_names.at(i),theta.at(i));
  }

  Rcpp::Rcout << " --- done initializing global variables and parameters --- " << std::endl;

  Rcpp::Rcout << " --- initializing simulation objects in memory --- " << std::endl;

  /* make the houses */
  for(size_t i=0; i<nhouse; i++){

    /* biting weight */
    double psi_i = Rcpp::as<double>(Rcpp::as<Rcpp::List>(house_param[i])["psi"]);

    // psi.emplace_back(0.);
    psi[i] = (double)psi_i;

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
    Rcpp::as<int>(mosy_param["K_eq"]),
    Rcpp::as<double>(mosy_param["lambda_v"])
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
        Rcpp::as<double>(hum_i["epsilon"]),
        Rcpp::as<double>(hum_i["lambda"]),
        Rcpp::as<double>(hum_i["phi"]),
        Rcpp::as<double>(hum_i["prDetectAMic"]),
        Rcpp::as<double>(hum_i["prDetectAPCR"]),
        Rcpp::as<double>(hum_i["prDetectUPCR"]),
        Rcpp::as<double>(hum_i["c"]),
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

    // run mosquito biting (mosy -> house transmission but to global)
    feeding_cycle(mosy_pop);

    // mosy -> house
    update_EIR(houses);

    // AFTER THIS POINT HUMANS/MOSY ARE CONDITIONALLY INDEPENDENT OF EACH OTHER

    // mosquito sim
    euler_step(mosy_pop);

    // human simulation functions
    one_day_update(houses);
    one_day_births(houses);
    one_day_deaths(houses);

    // bookkeeping before we move on
    pb.increment();
    tnow++;
  }

  Rcpp::Rcout << std::endl << " --- end simulation --- " << std::endl;

  houses.clear();

  // return output
  Rcpp::DataFrame state = Rcpp::DataFrame::create(
    Rcpp::Named("time") = Rcpp::wrap(time_out),
    Rcpp::Named("S") = Rcpp::wrap(state_S),
    Rcpp::Named("E") = Rcpp::wrap(state_E),
    Rcpp::Named("T") = Rcpp::wrap(state_T),
    Rcpp::Named("D") = Rcpp::wrap(state_D),
    Rcpp::Named("A") = Rcpp::wrap(state_A),
    Rcpp::Named("U") = Rcpp::wrap(state_U),
    Rcpp::Named("P") = Rcpp::wrap(state_P)
  );

  Rcpp::DataFrame age = Rcpp::DataFrame::create(
    Rcpp::Named("time") = Rcpp::wrap(time_out),
    Rcpp::Named("all") = Rcpp::wrap(num_All),
    Rcpp::Named("2_10") = Rcpp::wrap(num_2_10),
    Rcpp::Named("0_5") = Rcpp::wrap(num_0_5),
    Rcpp::Named("5_10") = Rcpp::wrap(num_5_10),
    Rcpp::Named("10_15") = Rcpp::wrap(num_10_15),
    Rcpp::Named("15+") = Rcpp::wrap(num_15Plus)
  );

  Rcpp::DataFrame clinic = Rcpp::DataFrame::create(
    Rcpp::Named("time") = Rcpp::wrap(time_out),
    Rcpp::Named("all") = Rcpp::wrap(cinc_All),
    Rcpp::Named("2_10") = Rcpp::wrap(cinc_2_10),
    Rcpp::Named("0_5") = Rcpp::wrap(cinc_0_5),
    Rcpp::Named("5_10") = Rcpp::wrap(cinc_5_10),
    Rcpp::Named("10_15") = Rcpp::wrap(cinc_10_15),
    Rcpp::Named("15+") = Rcpp::wrap(cinc_15Plus)
  );

  Rcpp::DataFrame mosy = Rcpp::DataFrame::create(
    Rcpp::Named("time") = Rcpp::wrap(time_out),
    Rcpp::Named("S") = Rcpp::wrap(mosy_S),
    Rcpp::Named("E") = Rcpp::wrap(mosy_E),
    Rcpp::Named("I") = Rcpp::wrap(mosy_I)
  );

  Rcpp::DataFrame trans = Rcpp::DataFrame::create(
    Rcpp::Named("time") = Rcpp::wrap(time_out),
    Rcpp::Named("lambda_v") = Rcpp::wrap(lambda_v)
  );

  return Rcpp::List::create(
    Rcpp::Named("state") = state,
    Rcpp::Named("age") = age,
    Rcpp::Named("clinical_incidence") = clinic,
    Rcpp::Named("mosy") = mosy,
    Rcpp::Named("trans") = trans
  );
};
