/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  April 2019
 #
 #  the house
*/

#include "house.hpp"

// other headers we need
#include "human.hpp"
#include "globals.hpp"

#include "stats.hpp"
#include "intervention.hpp"


/* ################################################################################
#   constructor & destructor
################################################################################ */

house::house(const size_t id_, stat_map* global_stat_, intervention_manager* int_mgr_) :
  id(id_), W(0.), Y(0.), Z(0.), C(0.), n(0), EIR(0), IRS(false), IRS_time_off(0.),
  global_stat(global_stat_), int_mgr(int_mgr_)
{};

house::~house(){};


/* ################################################################################
#   tracking output
################################################################################ */

void track_state(const house_vector& houses){

  time_out.at(tnow) = tnow;

  for(auto& hh : houses){
    for(auto& h : hh->humans){
      if(h->state.compare("S") == 0){
        state_S.at(tnow) += 1;
      } else if(h->state.compare("E") == 0){
        state_E.at(tnow) += 1;
      } else if(h->state.compare("T") == 0){
        state_T.at(tnow) += 1;
      } else if(h->state.compare("D") == 0){
        state_D.at(tnow) += 1;
      } else if(h->state.compare("A") == 0){
        state_A.at(tnow) += 1;
      } else if(h->state.compare("U") == 0){
        state_U.at(tnow) += 1;
      } else if(h->state.compare("P") == 0){
        state_P.at(tnow) += 1;
      } else {
        Rcpp::stop("incorrect state detected");
      }
    }
  }

};

void track_age(const house_vector& houses){

  for(auto& hh : houses){
    for(auto& h : hh->humans){

      if((h->age >= 2.) && (h->age < 10.)){
        num_2_10.at(tnow) += 1;
      }
      if(h->age < 5.) {
        num_0_5.at(tnow) += 1;
      } else if((h->age >= 5.) && (h->age < 10.)){
        num_5_10.at(tnow) += 1;
      } else if((h->age >= 10.) && (h->age < 15.)){
        num_10_15.at(tnow) += 1;
      } else if(h->age >= 15.){
        num_15Plus.at(tnow) += 1;
      }

      num_All.at(tnow) += 1;

    }
  }

};



/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

// for each house

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh){

  hh->C = 0.;

  for(auto& h : hh->humans){
    hh->C += hh->pi.at(h->id) * get_w(h) * h->c;
  }

};

// update W (net probability of successfuly feeding)
void update_W(house_ptr& hh){

  hh->W = 0.;

  for(auto& h : hh->humans){
    hh->W += hh->pi.at(h->id) * get_w(h);
  }

};

// update Z (net probability of being repelled w/out feeding)
void update_Z(house_ptr& hh){

  hh->Z = 0.;

  for(auto& h : hh->humans){
    hh->Z += hh->pi.at(h->id) * get_z(h);
  }

};

// normalize pi vector in a house
void normalize_pi(house_ptr& hh){

  const double pi_sum = std::accumulate(hh->pi.begin(), hh->pi.end(), 0.0,
    [](const double prev, const auto& element){
      return prev + element.second;
    }
  );

  std::for_each(hh->pi.begin(), hh->pi.end(), [pi_sum](auto& element){
    element.second /= pi_sum;
  });

};


// for the global landscape

// update global interface for mosquitos: THIS IS THE FUNCTION TO CALL (others are helpers)
// HUMAN -> MOSY

// this function updates (in order):
// pi for each house
// C for each house
// global CC
// global WW
// global ZZ
void update_biting(house_vector& houses){

  CC = 0.;
  WW = 0.;
  ZZ = 0.;

  for(auto& hh : houses){

    /* normalize biting weights (pi) */
    normalize_pi(hh);

    /* C for each house */
    update_C(hh);

    /* W for each house */
    update_W(hh);

    /* Z for each house */
    update_Z(hh);

    double psi_h = psi[hh->id];

    /* global CC is weighted average of household level C's */
    CC += psi_h * hh->C;

    /* global W is a weighted average of household level W's */
    WW += psi_h * hh->W;

    /* global Z is a weighted average of household level Z's */
    WW += psi_h * hh->Z;
  }

};

// MOSY -> HUMAN (after running feeding_cycle)
void update_EIR(house_vector& houses){

  for(size_t i=0; i<houses.size(); i++){
    houses[i]->EIR = EIR[i];
  }

};


/* ################################################################################
#   Interventions
################################################################################ */

void update_interventions_house(house_ptr& hh){
  if(hh->IRS && tnow > hh->IRS_time_off){
    hh->IRS = false;
  }
};

// spray the house
void apply_IRS(house_ptr& hh){

  double IRS_decay = parameters.at("IRS_decay");

  hh->IRS = true;
  hh->IRS_time_off = tnow + R::rgeom(IRS_decay);

};

// give ITNs to everyone in the house
void apply_ITN(house_ptr& hh){

  for(auto& h : hh->humans){
    give_ITN(h);
  }

};

// apply MDA; give drugs to everyone in the house
void apply_MDA(house_ptr& hh){

  for(auto& h : hh->humans){
    /* susceptibles, asymptomatic patent, and asymptomatic sub-patent go to P (chemoprophylaxis) */
    if(h->state.compare("S") == 0 || h->state.compare("A") == 0 || h->state.compare("U") == 0){
      h->state = "P";
    /* incubating and untreated clinical episodes go to T (treated clinical) */
    } else if(h->state.compare("E") == 0 || h->state.compare("D") == 0){
      h->state = "T";
    }
  }

};

// apply RACD; test people via PCR, only give drugs to those who test positive
void apply_RACD_PCR(house_ptr& hh){

  for(auto& h : hh->humans){

    /* untreated clinical cases are always caught */
    if(h->state.compare("D") == 0){
      h->state = "T";
    /* patent asymptomatic cases: roll the dice */
    } else if(h->state.compare("A") == 0){
      /* PCR test on the person */
      if(R::runif(0.,1.) < h->prDetectAPCR){
        h->state = "P";
      }
    /* sub-patent asymptomatic cases: roll the dice (just for PCR or LAMP) */
    } else if(h->state.compare("U") == 0){
      /* PCR test on the person */
      if(R::runif(0.,1.) < h->prDetectUPCR){
        h->state = "P";
      }
    }

  }

};

// apply RACD; test people via microscopy, only give drugs to those who test positive
void apply_RACD_Mic(house_ptr& hh){

  for(auto& h : hh->humans){

    /* untreated clinical cases are always caught */
    if(h->state.compare("D") == 0){
      h->state = "T";
    /* patent asymptomatic cases: roll the dice */
    } else if(h->state.compare("A") == 0){
      /* PCR test on the person */
      if(R::runif(0.,1.) < h->prDetectAMic){
        h->state = "P";
      }
    }
  }

  /* sub-patent asymptomatic cases are invisible to microscopy */

};

// apply RACD; test people via LAMP, only give drugs to those who test positive
void apply_RACD_LAMP(house_ptr& hh){
  Rcpp::stop("can't use RACD with LAMP for Pf detection yet!");
};


/* ################################################################################
#   demographics
################################################################################ */

/* shamelessly "referenced" from Knuth, The Art of Computer Programming Vol 2, section 4.2.2 */
// the mean immunity among 18-22 year olds
double mean_ICA18_22(house_vector& houses){
  double avg = 0.;
  int t = 1;
  for (auto& hh : houses) {
    for(auto& h : hh->humans){
      if((h->age >= 18.) && (h->age < 22.)){
        avg += (h->ICA - avg) / t;
        ++t;
      }
    }
  }
  // weird hack if there aren't any 18-22 yr olds
  // just take avg ICA of entire pop
  if(avg < 0.001){
    avg = 0.;
    t = 1;
    for (auto& hh : houses) {
      for(auto& h : hh->humans){
        avg += (h->ICA - avg) / t;
        ++t;
      }
    }
  }
  return avg;
};

// bring out yer dead!
void one_day_deaths(house_vector& houses){
  for(auto& hh : houses){
    hh->humans.remove_if(
      [](auto& h) {return !h->alive;}
    );
  }
};

// the respawn point
void one_day_births(house_vector& houses){

  size_t hpop = num_All.at(tnow);
  double mu = parameters.at("mu");

  size_t nbirth = (size_t)R::rbinom((double)hpop, mu);

  if(nbirth > 0){

    double ICA18_22 = mean_ICA18_22(houses);

    double sigma2 = parameters.at("sigma2");

    double PM = parameters.at("PM");

    double phi0 = parameters.at("phi0");
    double phi1 = parameters.at("phi1");
    double kappaC = parameters.at("kappaC");
    double IC0 = parameters.at("IC0");

    for(size_t i=0; i<nbirth; i++){

      // put newborns in the smallest houses for ... reasons
      size_t smallest_hh = std::distance(houses.begin(),std::min_element(houses.begin(), houses.end(), [](auto& hh1, auto& hh2) {
        return hh1->n < hh2->n;
      }));

      // sample this person's biting heterogeneity
      double zeta = R::rlnorm(-sigma2/2., std::sqrt(sigma2));

      // P(clinical disease | infection)
      double phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(PM*ICA18_22/IC0,kappaC))));

      // put the human in their new home
      houses.at(smallest_hh)->humans.emplace_back(std::make_unique<human>(
        0.,
        houses.at(smallest_hh).get(),
        zeta,
        0.,
        0.,
        0.,
        (PM * ICA18_22),
        0.,
        0.,
        phi,
        1.,
        1.,
        1.,
        0.,
        "S"
      ));

    }

  }

};


/* ################################################################################
#   daily simulation
################################################################################ */

// update the dynamics of a single dwelling
void one_day_update_house(house_ptr& hh){

  for(auto& h : hh->humans){
    one_day_update_human(h);
  }

};

// the update for all humans/dwellings
void one_day_update(house_vector& houses){

  for(auto& hh : houses){
    // simulate the humans
    one_day_update_house(hh);
    // update interventions at this house
    update_interventions_house(hh);
  }

};
