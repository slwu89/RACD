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
 #  the human
*/

#include "human.hpp"

// other headers we need
#include "house.hpp"
#include "globals.hpp"
#include "stats.hpp"
#include "intervention.hpp"

// for ageing
static const double one_day = 1./365.;


/* ################################################################################
#   constructor & destructor
################################################################################ */

human::human(const double age_,
      house* house_ptr_,
      const double zeta_,
      const double IB_,
      const double ID_,
      const double ICA_,
      const double ICM_,
      const double epsilon_,
      const double lambda_,
      const double phi_,
      const double prDetectAMic_,
      const double prDetectAPCR_,
      const double prDetectUPCR_,
      const double c_,
      const std::string state_) :
      id(globals::instance().get_global_hid()),
      age(age_),
      alive(true),
      house_ptr(house_ptr_),
      zeta(zeta_),
      IB(IB_),
      ID(ID_),
      ICA(ICA_),
      ICM(ICM_),
      epsilon(epsilon_),
      lambda(lambda_),
      phi(phi_),
      prDetectAMic(prDetectAMic_),
      prDetectAPCR(prDetectAPCR_),
      prDetectUPCR(prDetectUPCR_),
      c(c_),
      state(state_),
      days_latent(0),
      ITN(false),
      ITN_time_off(0.)
{
  // add my biting to the hash table
  double a0 = globals::instance().get_pmap().at("a0");
  double rho = globals::instance().get_pmap().at("rho");
  double pi = zeta * (1. - rho * std::exp(-age/a0));
  house_ptr->pi.emplace(id,pi);

  // let the house know i showed up
  house_ptr->n += 1;

};

human::~human(){};


/* ################################################################################
#   individual level tracking
################################################################################ */

// track clinical incidence
// this function does double duty as it also lets the
// intervention manager know a case popped up
void track_cinc(const human_ptr& h){

  size_t j;

  if((h->age >= 2.) && (h->age < 10.)){
    globals::instance().push_cinc_2_10();
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

  globals::instance().push_cinc_age(j);

  // let the intervention mgr know about this person's case
  h->house_ptr->int_mgr->add_cinc(h->house_ptr->id);

};


/* ################################################################################
#   State transitions for our little Markov humans
################################################################################ */

void mortality(human_ptr& human){
  double randNum = R::runif(0.0,1.0);
  double mu = globals::instance().get_pmap().at("mu");
  if(randNum <= mu){
    human->alive = false;
    human->house_ptr->n -= 1;
    remove_pi(human);
  }
}

/* S: susceptible */
void S_compartment(human_ptr& human){

  double randNum = R::runif(0.0,1.0);

  if(randNum <= human->lambda){
    human->state = "E";
    human->days_latent = 0;
  }
};


/* E: latent period */
void E_compartment(human_ptr& human){

  double dE = globals::instance().get_pmap().at("dE");
  double fT = globals::instance().get_pmap().at("fT");

  if(human->days_latent < dE){
    human->days_latent++;
  } else {
    double randNum = R::runif(0.0,1.0);

    double phi = human->phi;

    // Treated clinical infection (E -> T)
    if(randNum <= phi*fT){
      human->state = "T";
      human->days_latent = 0;

      track_cinc(human);
    }

    // Untreated clinical infection (E -> D)
    if((randNum > phi*fT) && (randNum <= phi)){
      human->state = "D";
      human->days_latent = 0;

      track_cinc(human);
    }

    // Asymptomatic infection (E -> A)
    if(randNum > phi){
      human->state = "A";
      human->days_latent = 0;
    }
  }
};

/* T: treated clinical disease */
void T_compartment(human_ptr& human){

  double dT = globals::instance().get_pmap().at("dT");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dT)){
    human->state = "P";
  }

};

/* D: untreated clinical disease */
void D_compartment(human_ptr& human){

  double dD = globals::instance().get_pmap().at("dD");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dD)){
    human->state = "A";
  }
};

/* A: asymptomatic patent (detectable by microscopy) infection */
void A_compartment(human_ptr& human){

  double fT = globals::instance().get_pmap().at("fT");
  double dA = globals::instance().get_pmap().at("dA");

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (A -> T)
  if(randNum <= phi*fT*lambda){
    human->state = "T";

    track_cinc(human);
  }
  // Untreated clinical infection (A -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";

    track_cinc(human);
  }
  // Progression to asymptomatic sub-patent infection (A -> U):
  if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1.0/dA)))) {
    human->state = "U";
  }
};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void U_compartment(human_ptr& human){

  double fT = globals::instance().get_pmap().at("fT");
  double dU = globals::instance().get_pmap().at("dU");

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (U -> T):
  if(randNum <= phi*fT*lambda){
    human->state = "T";

    track_cinc(human);
  }

  // Untreated clinical infection (U -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";

    track_cinc(human);
  }

  // Asymptomatic infection (U -> A)
  if((randNum > phi*lambda) && (randNum <= lambda)){
    human->state = "A";
  }

  // Recovery to susceptible (U -> S):
  if((randNum > lambda) && (randNum <= (lambda + (1.0/dU)))){
    human->state = "S";
  }
};

/* P: protection due to chemoprophylaxis treatment */
void P_compartment(human_ptr& human){

  double dP = globals::instance().get_pmap().at("dP");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dP)){
    human->state = "S";
  }
};


/* ################################################################################
#   Immunity Functions
################################################################################ */

/* immunity */
/* Update values for:
 * 1. Pre-erythrocytic immunity (IB, reduces the probability of infection
 *  following an infectious challenge)
 * 2. Acquired clinical immunity (ICA, reduces the probability of clinical
 *  disease, acquired from previous exposure)
 * 3. Maternal clinical immunity (ICM, reduces the probability of clinical
 *  disease, acquired maternally)
 * 4. Detection immunity (ID, a.k.a. blood-stage immunity, reduces the
 *  probability of detection and reduces infectiousness to mosquitoes)
 */
void update_immunity(human_ptr& human){

  double uB = globals::instance().get_pmap().at("uB");
  double uC = globals::instance().get_pmap().at("uC");
  double uD = globals::instance().get_pmap().at("uD");
  double dB = globals::instance().get_pmap().at("dB");
  double dC = globals::instance().get_pmap().at("dC");
  double dID = globals::instance().get_pmap().at("dID");
  double dM = globals::instance().get_pmap().at("dM");

  double epsilon = human->epsilon;
  double lambda = human->lambda;

  double IB = human->IB;
  double ICA = human->ICA;
  double ICM = human->ICM;
  double ID = human->ID;

  human->IB = IB + (epsilon/(epsilon*uB + 1.)) - (IB)*(1./dB);
  human->ICA = ICA + (lambda/(lambda*uC + 1.)) - (ICA)*(1./dC);
  human->ICM = ICM - (ICM)*(1./dM);
  human->ID = ID + (lambda/(lambda*uD + 1.)) - (ID)*(1./dID);

};

/* lambda */
// psi is the psi of my house
void update_lambda(human_ptr& human){

  double b0 = globals::instance().get_pmap().at("b0");
  double b1 = globals::instance().get_pmap().at("b1");
  double IB0 = globals::instance().get_pmap().at("IB0");
  double kappaB = globals::instance().get_pmap().at("kappaB");
  double IB = human->IB;

  // my EIR (house EIR * P(it bites me))
  // double EIR_h = EIR.at(human->house_ptr->id) * human->house_ptr->pi.at(human->id);
  double EIR_h = globals::instance().get_EIR().at(human->house_ptr->id) * human->house_ptr->pi.at(human->id);

  human->epsilon = EIR_h * get_y(human); // term to account for possible effect of intervention
  double b = b0*(b1 + ((1.-b1)/(1. + std::pow((IB/IB0),kappaB))));
  human->lambda = EIR_h * b;

  // track EIR values
  globals::instance().get_stats().at("EIR").Push(human->epsilon);

  // track b values
  globals::instance().get_stats().at("b").Push(b);
};

/* phi */
void update_phi(human_ptr& human){

  double phi0 = globals::instance().get_pmap().at("phi0");
  double phi1 = globals::instance().get_pmap().at("phi1");
  double IC0 = globals::instance().get_pmap().at("IC0");
  double kappaC = globals::instance().get_pmap().at("kappaC");

  double ICA = human->ICA;
  double ICM = human->ICM;

  human->phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
void update_q(human_ptr& human){

  double fD0 = globals::instance().get_pmap().at("fD0");
  double aD = globals::instance().get_pmap().at("aD");
  double gammaD = globals::instance().get_pmap().at("gammaD");
  double d1 = globals::instance().get_pmap().at("d1");
  double ID0 = globals::instance().get_pmap().at("ID0");
  double kappaD = globals::instance().get_pmap().at("kappaD");
  double alphaA = globals::instance().get_pmap().at("alphaA");
  double alphaU = globals::instance().get_pmap().at("alphaU");

  double ID = human->ID;
  double age = human->age;

  double fD = 1. - ((1. - fD0)/(1. + std::pow((age/aD),gammaD)));
  double q = d1 + ((1. - d1)/(1. + (fD*std::pow((ID/ID0),kappaD))*fD));

  human->prDetectAMic = q;
  human->prDetectAPCR = std::pow(q,alphaA);
  human->prDetectUPCR = std::pow(q,alphaU);

};


/* ################################################################################
#   Infectiousness to Mosquitoes
################################################################################ */

void infectiousness_S(human_ptr& human){
  human->c = 0.;
};

void infectiousness_E(human_ptr& human){
  human->c = 0.;
};

void infectiousness_T(human_ptr& human){
  human->c = globals::instance().get_pmap().at("cT");
};

void infectiousness_D(human_ptr& human){
  human->c = globals::instance().get_pmap().at("cD");
};

void infectiousness_A(human_ptr& human){
  double cU = globals::instance().get_pmap().at("cU");
  double cD = globals::instance().get_pmap().at("cD");
  double gammaI = globals::instance().get_pmap().at("gammaI");
  human->c = cU + (cD - cU)*std::pow(human->prDetectAMic,gammaI);
};

void infectiousness_U(human_ptr& human){
  human->c = globals::instance().get_pmap().at("cU");
};

void infectiousness_P(human_ptr& human){
  human->c = 0.;
};


/* ################################################################################
#   Mosquito Approch probabilities (what happens when a bloodsucker tries to bite me)
################################################################################ */

/* probability of successful biting */
double get_w(human_ptr& human){

  bool IRS = human->house_ptr->IRS;
  bool ITN = human->ITN;

  /* none */
  if(!ITN && !IRS){
    return 1.0;
  /* IRS only */
  } else if(IRS && !ITN){

    double phiI = globals::instance().get_pmap().at("phiI");
    double rS = globals::instance().get_pmap().at("rIRS");
    double sS = globals::instance().get_pmap().at("sIRS");

    return (1. - phiI) + (phiI * (1. - rS) * sS);
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = globals::instance().get_pmap().at("phiB");
    double sN = globals::instance().get_pmap().at("sITN");

    return (1. - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = globals::instance().get_pmap().at("phiI");
    double rS = globals::instance().get_pmap().at("rIRS");
    double sS = globals::instance().get_pmap().at("sIRS");

    double phiB = globals::instance().get_pmap().at("phiB");
    double sN = globals::instance().get_pmap().at("sITN");

    return (1. - phiI) + (phiB * (1. - rS) * sN * sS) + ((phiI - phiB) * (1. - rS) * sS);
  } else {
    Rcpp::stop("error: invalid combination of ITN/IRS");
  }
};

/* probability of biting */
double get_y(human_ptr& human){

  bool IRS = human->house_ptr->IRS;
  bool ITN = human->ITN;

  /* none */
  if(!ITN && !IRS){
    return 1.0;
  /* IRS only */
  } else if(IRS && !ITN){

    double phiI = globals::instance().get_pmap().at("phiI");
    double rS = globals::instance().get_pmap().at("rIRS");

    return (1.0 - phiI) + (phiI * (1. - rS));
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = globals::instance().get_pmap().at("phiB");
    double sN = globals::instance().get_pmap().at("sITN");

    return (1.0 - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = globals::instance().get_pmap().at("phiI");
    double rS = globals::instance().get_pmap().at("rIRS");

    double phiB = globals::instance().get_pmap().at("phiB");
    double sN = globals::instance().get_pmap().at("sITN");

    return (1.0 - phiI) + (phiB * (1. - rS) * sN) + ((phiI - phiB) * (1. - rS) );
  } else {
    Rcpp::stop("error: invalid combination of ITN/IRS");
  }
};

/* probability of repellency*/
double get_z(human_ptr& human){

  bool IRS = human->house_ptr->IRS;
  bool ITN = human->ITN;

  /* none */
  if(!ITN && !IRS){
    return 0.0;
  /* IRS only */
  } else if(IRS && !ITN){

    double phiI = globals::instance().get_pmap().at("phiI");
    double rS = globals::instance().get_pmap().at("rIRS");

    return phiI * rS;
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = globals::instance().get_pmap().at("phiB");
    double rN = globals::instance().get_pmap().at("rN");

    return phiB * rN;
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = globals::instance().get_pmap().at("phiI");
    double rS = globals::instance().get_pmap().at("rIRS");

    double phiB = globals::instance().get_pmap().at("phiB");
    double rN = globals::instance().get_pmap().at("rN");

    return (phiB * (1. - rS) * rN) + (phiI * rS);
  } else {
    Rcpp::stop("error: invalid combination of ITN/IRS");
  }
};


/* ################################################################################
#   bookkeeping (bites and biting weights)
#   these work on humans and houses
################################################################################ */

// add my biting weight to the hash table
void add_pi(human_ptr& human){

  double a0 = globals::instance().get_pmap().at("a0");
  double rho = globals::instance().get_pmap().at("rho");
  double pi = human->zeta * (1. - rho * std::exp(-human->age/a0));

  // put it into the hash table
  human->house_ptr->pi.emplace(human->id,pi);

};

// take out my biting weight
void remove_pi(human_ptr& human){
  human->house_ptr->pi.erase(human->id);
}

// update pi (do this daily because I age)
void update_pi(human_ptr& human){

  double a0 = globals::instance().get_pmap().at("a0");
  double rho = globals::instance().get_pmap().at("rho");
  double pi = human->zeta * (1. - rho * std::exp(-human->age/a0));

  // update the hash table
  human->house_ptr->pi.at(human->id) = pi;

};


/* ################################################################################
#   Humans: daily update
################################################################################ */

void one_day_update_human(human_ptr& human){

  // mortality
  mortality(human);

  // if they survived the call of the beyond
  if(human->alive){

    // state update
    state_functions.at(human->state)(human);

    // update infectiousness to mosquitos
    infectivity_functions.at(human->state)(human);

    // update age
    human->age += one_day;

    // update immunity
    update_immunity(human);

    // update lambda
    update_lambda(human);

    // update phi
    update_phi(human);

    // update q
    update_q(human);

    // update pi
    update_pi(human);

    // update interventions
    update_interventions_human(human);

  }

};


/* ################################################################################
#   interventions
################################################################################ */

// called before exiting daily update; check if interventions expire
void update_interventions_human(human_ptr& human){

  if(human->ITN && globals::instance().get_tnow() >= human->ITN_time_off){
    human->ITN = false;
  }

};

void give_ITN(human_ptr& human){

  double ITN_decay = globals::instance().get_pmap().at("ITN_decay");

  human->ITN = true;
  human->ITN_time_off = globals::instance().get_tnow() + R::rgeom(ITN_decay);

};
