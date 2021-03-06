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

house::house(const size_t id_, intervention_manager* int_mgr_) :
  id(id_), W(0.), Y(0.), Z(0.), C(0.), n(0), EIR(0), IRS(false), IRS_time_off(0.),
  int_mgr(int_mgr_)
{};

house::~house(){};


/* ################################################################################
#   tracking output
################################################################################ */

void track_state_age(const house_vector& houses){

  size_t i;
  size_t j;

  for(auto& hh : houses){
    for(auto& h : hh->humans){

      // get their state
      if(h->state.compare("S") == 0){
        i = 0;
      } else if(h->state.compare("E") == 0){
        i = 1;
      } else if(h->state.compare("T") == 0){
        i = 2;
      } else if(h->state.compare("D") == 0){
        i = 3;
      } else if(h->state.compare("A") == 0){
        i = 4;
      } else if(h->state.compare("U") == 0){
        i = 5;
      } else if(h->state.compare("P") == 0){
        i = 6;
      } else {
        Rcpp::stop("incorrect state detected");
      }

      // get their age category
      if((h->age >= 2.) && (h->age < 10.)){
        globals::instance().push_state_2_10_tnow(i);
      }

      if(h->age < 5.) {
        j = 2;
      } else if((h->age >= 5.) && (h->age < 10.)){
        j = 3;
      } else if((h->age >= 10.) && (h->age < 15.)){
        j = 4;
      } else if(h->age >= 15.){
        j = 5;
      } else {
        Rcpp::stop("invalid age for human");
      }

      // log output (once for all, and then for my age-class)
      globals::instance().push_state_age_tnow(i,0);
      globals::instance().push_state_age_tnow(i,j);
    }
  }


};



/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

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

  globals::instance().zero_CC();
  globals::instance().zero_WW();
  globals::instance().zero_ZZ();

  for(auto& hh : houses){

    /* normalize biting weights (pi) */
    normalize_pi(hh);

    /* W for each house */
    update_W(hh);

    /* C for each house */
    update_C(hh);

    /* Z for each house */
    update_Z(hh);

    double psi_h = globals::instance().get_psi().at(hh->id);

    /* global CC is weighted average of household level C's */
    globals::instance().inc_CC(psi_h * hh->C);

    /* global W is a weighted average of household level W's */
    globals::instance().inc_WW(psi_h * hh->W);

    /* global Z is a weighted average of household level Z's */
    globals::instance().inc_ZZ(psi_h * hh->Z);

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

// update W (net probability of successfuly feeding)
void update_W(house_ptr& hh){

  hh->W = 0.;

  for(auto& h : hh->humans){
    hh->W += hh->pi.at(h->id) * get_w(h);
  }

};

// for each house

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh){

  hh->C = 0.;

  for(auto& h : hh->humans){
    hh->C += h->c * (hh->pi.at(h->id) * get_w(h)) / hh->W;
  }

};

// update Z (net probability of being repelled w/out feeding)
void update_Z(house_ptr& hh){

  hh->Z = 0.;

  for(auto& h : hh->humans){
    hh->Z += hh->pi.at(h->id) * get_z(h);
  }

};

// MOSY -> HUMAN (after running feeding_cycle)
void update_EIR(house_vector& houses){

  for(size_t i=0; i<houses.size(); i++){
    // houses[i]->EIR = EIR[i];
    houses[i]->EIR = globals::instance().get_EIR().at(i);
  }

};


/* ################################################################################
#   Interventions
################################################################################ */

void update_interventions_house(house_ptr& hh){
  if(hh->IRS && globals::instance().get_tnow() > hh->IRS_time_off){
    hh->IRS = false;
  }
};

// spray the house
void apply_IRS(house_ptr& hh){

  double IRS_decay = globals::instance().get_pmap().at("IRS_decay");

  hh->IRS = true;
  hh->IRS_time_off = globals::instance().get_tnow() + R::rgeom(IRS_decay);

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

  // BIRTHS
  int hpop = std::accumulate( globals::instance().get_state_age_tnow().column(0).begin(),  globals::instance().get_state_age_tnow().column(0).end(),0);
  // size_t hpop = num_All.at(tnow);
  double mu = globals::instance().get_pmap().at("mu");

  size_t nbirth = (size_t)R::rbinom((double)hpop, mu);

  if(nbirth > 0){

    double ICA18_22 = mean_ICA18_22(houses);

    double sigma2 = globals::instance().get_pmap().at("sigma2");

    double PM = globals::instance().get_pmap().at("PM");

    double phi0 = globals::instance().get_pmap().at("phi0");
    double phi1 = globals::instance().get_pmap().at("phi1");
    double kappaC = globals::instance().get_pmap().at("kappaC");
    double IC0 = globals::instance().get_pmap().at("IC0");

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

  // IMPORTATION
  double import_rate = globals::instance().get_pmap().at("import_rate");
  int import_case = (int)R::rpois(import_rate);

  if(import_case > 0){

    // grab the data we need to make humans
    double ICA18_22 = mean_ICA18_22(houses);
    double sigma2 = globals::instance().get_pmap().at("sigma2");
    double PM = globals::instance().get_pmap().at("PM");
    double phi0 = globals::instance().get_pmap().at("phi0");
    double phi1 = globals::instance().get_pmap().at("phi1");
    double kappaC = globals::instance().get_pmap().at("kappaC");
    double IC0 = globals::instance().get_pmap().at("IC0");

    // fixed params
    double IB_imp = globals::instance().get_pmap().at("IB_imp");
    double ID_imp = globals::instance().get_pmap().at("ID_imp");
    double ICA_imp = globals::instance().get_pmap().at("ICA_imp");

    // put each person in a house
    for(int i=0; i<import_case; i++){

      // randomly sample destination
      int dest = globals::instance().sample_lookup();

      // sample this person's biting heterogeneity
      double zeta = R::rlnorm(-sigma2/2., std::sqrt(sigma2));

      // P(clinical disease | infection)
      double phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(PM*ICA18_22/IC0,kappaC))));

      // put the infected human in their new home
      houses.at(dest)->humans.emplace_back(std::make_unique<human>(
        21,
        houses.at(dest).get(),
        zeta,
        IB_imp,
        ID_imp,
        ICA_imp,
        (PM * ICA18_22),
        0.,
        0.,
        phi,
        1.,
        1.,
        1.,
        0.,
        "D"
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
