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
      id(global_hid),
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
  // after we take our id, increment for the next person!
  global_hid++;

  // add my biting to the hash table
  double a0 = parameters.at("a0");
  double rho = parameters.at("rho");
  double pi = zeta * (1. - rho * std::exp(-age/a0));
  house_ptr->pi.emplace(id,pi);

  // let the house know i showed up
  house_ptr->n += 1;

};

human::~human(){
  // Rcpp::Rcout << "human: " << id << " is dying because they are dead: " << alive << " in house: " << house_ptr->id << "\n";
};


/* ################################################################################
#   individual level tracking
################################################################################ */

// track clinical incidence
void track_cinc(const human_ptr& h){

  cinc_All.at(tnow) += 1;

  if((h->age >= 2.) && (h->age < 10.)){
    cinc_2_10.at(tnow) += 1;
  }
  if(h->age < 5.) {
    cinc_0_5.at(tnow) += 1;
  } else if((h->age >= 5.) && (h->age < 10.)){
    cinc_5_10.at(tnow) += 1;
  } else if((h->age >= 10.) && (h->age < 15.)){
    cinc_10_15.at(tnow) += 1;
  } else if(h->age >= 15.){
    cinc_15Plus.at(tnow) += 1;
  }

};


/* ################################################################################
#   State transitions for our little Markov humans
################################################################################ */

void mortality(human_ptr& human){
  double randNum = R::runif(0.0,1.0);
  double mu = parameters.at("mu");
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

  double dE = parameters.at("dE");
  double fT = parameters.at("fT");

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

  double dT = parameters.at("dT");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dT)){
    human->state = "P";
  }

};

/* D: untreated clinical disease */
void D_compartment(human_ptr& human){

  double dD = parameters.at("dD");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dD)){
    human->state = "A";
  }
};

/* A: asymptomatic patent (detectable by microscopy) infection */
void A_compartment(human_ptr& human){

  double fT = parameters.at("fT");
  double dA = parameters.at("dA");

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

  double fT = parameters.at("fT");
  double dU = parameters.at("dU");

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

  double dP = parameters.at("dP");
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

  double uB = parameters.at("uB");
  double uC = parameters.at("uC");
  double uD = parameters.at("uD");
  double dB = parameters.at("dB");
  double dC = parameters.at("dC");
  double dID = parameters.at("dID");
  double dM = parameters.at("dM");

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

  double b0 = parameters.at("b0");
  double b1 = parameters.at("b1");
  double IB0 = parameters.at("IB0");
  double kappaB = parameters.at("kappaB");
  double IB = human->IB;

  // my EIR (house EIR * P(it bites me))
  double EIR_h = EIR.at(human->house_ptr->id) * human->house_ptr->pi.at(human->id);

  human->epsilon = EIR_h * get_y(human); // term to account for possible effect of intervention
  double b = b0*(b1 + ((1.-b1)/(1. + std::pow((IB/IB0),kappaB))));
  human->lambda = EIR_h * b;

  // track EIR values
  human->house_ptr->global_stat->at("EIR")->Push(human->epsilon);

};

/* phi */
void update_phi(human_ptr& human){

  double phi0 = parameters.at("phi0");
  double phi1 = parameters.at("phi1");
  double IC0 = parameters.at("IC0");
  double kappaC = parameters.at("kappaC");

  double ICA = human->ICA;
  double ICM = human->ICM;

  human->phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
void update_q(human_ptr& human){

  double fD0 = parameters.at("fD0");
  double aD = parameters.at("aD");
  double gammaD = parameters.at("gammaD");
  double d1 = parameters.at("d1");
  double ID0 = parameters.at("ID0");
  double kappaD = parameters.at("kappaD");
  double alphaA = parameters.at("alphaA");
  double alphaU = parameters.at("alphaU");

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
  human->c = parameters.at("cT");
};

void infectiousness_D(human_ptr& human){
  human->c = parameters.at("cD");
};

void infectiousness_A(human_ptr& human){
  double cU = parameters.at("cU");
  double cD = parameters.at("cD");
  double gammaI = parameters.at("gammaI");
  human->c = cU + (cD - cU)*std::pow(human->prDetectAMic,gammaI);
};

void infectiousness_U(human_ptr& human){
  human->c = parameters.at("cU");
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

    double phiI = parameters.at("phiI");
    double rS = parameters.at("rIRS");
    double sS = parameters.at("sIRS");

    return (1. - phiI) + (phiI * (1. - rS) * sS);
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = parameters.at("phiB");
    double sN = parameters.at("sITN");

    return (1. - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = parameters.at("phiI");
    double rS = parameters.at("rIRS");
    double sS = parameters.at("sIRS");

    double phiB = parameters.at("phiB");
    double sN = parameters.at("sITN");

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

    double phiI = parameters.at("phiI");
    double rS = parameters.at("rIRS");

    return (1.0 - phiI) + (phiI * (1. - rS));
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = parameters.at("phiB");
    double sN = parameters.at("sITN");

    return (1.0 - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = parameters.at("phiI");
    double rS = parameters.at("rIRS");

    double phiB = parameters.at("phiB");
    double sN = parameters.at("sITN");

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

    double phiI = parameters.at("phiI");
    double rS = parameters.at("rIRS");

    return phiI * rS;
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = parameters.at("phiB");
    double rN = parameters.at("rN");

    return phiB * rN;
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = parameters.at("phiI");
    double rS = parameters.at("rIRS");

    double phiB = parameters.at("phiB");
    double rN = parameters.at("rN");

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

  double a0 = parameters.at("a0");
  double rho = parameters.at("rho");
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

  double a0 = parameters.at("a0");
  double rho = parameters.at("rho");
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
    human->age += (1./365.);

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

  if(human->ITN && tnow >= human->ITN_time_off){
    human->ITN = false;
  }

};

void give_ITN(human_ptr& human){

  double ITN_decay = parameters.at("ITN_decay");

  human->ITN = true;
  human->ITN_time_off = R::rgeom(ITN_decay);

};
