/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  November 2019
 #
 #  Mosquitoes
*/

#include "house.hpp"
#include "human.hpp"

// [[Rcpp::plugins(cpp14)]]


/* ################################################################################
#   house ctor/dtor
################################################################################ */

house::house() :
  id(global_id++), W(0.), Y(0.), Z(0.), C(0.), n(0), EIR(0), IRS(false), IRS_deploy(0), IRS_decay(0), cinc(0)
{
  Rcpp::Rcout << "house id: " << id << ", ctor called at " << this << "\n";
};

house::~house(){
  Rcpp::Rcout << "house id: " << id << ", dtor called at " << this << "\n";
};


/* ################################################################################
#   initialize houses and return to R
################################################################################ */

Rcpp::XPtr<house_vector> init_houses(
  const Rcpp::List& humans_param,
  const Rcpp::List& house_param
){

};


// cinc needs to be reset daily before anything starts

/* ################################################################################
#   Biting stuff (the interface btwn humans and mosquitos)
################################################################################ */

// update global interface for mosquitos: THIS IS THE FUNCTION TO CALL (others are helpers)
// HUMAN -> MOSY
// returns bite_probs: vector of WW,ZZ,CC for mosquitoes
std::vector<double> update_biting(house_vector& houses){

  std::vector<double> biting(3,0.); // WW,ZZ,CC

  for(auto& hh : houses){

    /* normalize biting weights (pi) */
    normalize_pi(hh);

    /* update W */
    update_W(hh);

    /* update C (kappa) */
    update_C(hh);

    /* update Z */
    update_Z(hh);

    /* update output */
    double psi = hh->psi;

    biting[0] += psi * hh->W;
    biting[1] += psi * hh->C;
    biting[2] += psi * hh->Z;

  }

  return biting;
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

// update C (net infectivity of this house to mosquitos)
void update_C(house_ptr& hh){

  hh->C = 0.;

  for(auto& h: hh->humans){
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
// EIR is the vector from mosquitoes
void update_EIR(house_vector& houses, const std::vector<double>& EIR){
  int i = 0;
  for(auto& h : houses){
    h->EIR = EIR[i];
    h->cinc = 0;
    i++;
  }
};


/* ################################################################################
#   Interventions
################################################################################ */

void update_interventions_house(house_ptr& hh, const int tnow){
  if(hh->IRS && tnow > hh->IRS_decay){
    hh->IRS = false;
  }
};

// spray the house
void apply_IRS(house_ptr& hh, const int tnow){

  double IRS_decay = hh->par_ptr->at("IRS_decay");

  hh->IRS = true;
  hh->IRS_deploy = tnow;
  hh->IRS_decay = tnow + R::rgeom(IRS_decay);

};

// give ITNs to everyone in the house
void apply_ITN(house_ptr& hh, const int tnow){

  for(auto& h : hh->humans){
    give_ITN(h,tnow);
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
