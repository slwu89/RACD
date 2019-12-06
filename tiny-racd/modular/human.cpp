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

#include "house.hpp"

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
      const std::string state_
    ) :
      id(global_id++),
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
      ITN_given(0),
      ITN_decay(0)
{
  double rho = house_ptr->par_ptr->at("rho");
  double a0 = house_ptr->par_ptr->at("a0");

  // add my biting to the hash table
  double pi = zeta * (1. - rho * std::exp(-age/a0));
  house_ptr->pi.emplace(id,pi);

  // let the house know i showed up
  house_ptr->n += 1;

};

human::~human(){};


/* ################################################################################
#   State transitions for our little Markov humans
################################################################################ */

void mortality(human_ptr& human){
  double randNum = R::runif(0.0,1.0);
  double mu = human->house_ptr->par_ptr->at("mu");
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

  double dE = human->house_ptr->par_ptr->at("dE");
  double fT = human->house_ptr->par_ptr->at("fT");

  if(human->days_latent < dE){
    human->days_latent++;
  } else {
    double randNum = R::runif(0.0,1.0);

    double phi = human->phi;

    // Treated clinical infection (E -> T)
    if(randNum <= phi*fT){
      human->state = "T";
      human->days_latent = 0;

      human->house_ptr->cinc += 1;
    }

    // Untreated clinical infection (E -> D)
    if((randNum > phi*fT) && (randNum <= phi)){
      human->state = "D";
      human->days_latent = 0;

      human->house_ptr->cinc += 1;
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

  double dT = human->house_ptr->par_ptr->at("dT");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dT)){
    human->state = "P";
  }

};

/* D: untreated clinical disease */
void D_compartment(human_ptr& human){

  double dD = human->house_ptr->par_ptr->at("dD");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dD)){
    human->state = "A";
  }
};

/* A: asymptomatic patent (detectable by microscopy) infection */
void A_compartment(human_ptr& human){

  double fT = human->house_ptr->par_ptr->at("fT");
  double dA = human->house_ptr->par_ptr->at("dA");

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (A -> T)
  if(randNum <= phi*fT*lambda){
    human->state = "T";
    human->house_ptr->cinc += 1;
  }
  // Untreated clinical infection (A -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";
    human->house_ptr->cinc += 1;
  }
  // Progression to asymptomatic sub-patent infection (A -> U):
  if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1.0/dA)))) {
    human->state = "U";
  }
};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void U_compartment(human_ptr& human){

  double fT = human->house_ptr->par_ptr->at("fT");
  double dU = human->house_ptr->par_ptr->at("dU");

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (U -> T):
  if(randNum <= phi*fT*lambda){
    human->state = "T";

    human->house_ptr->cinc += 1;
  }

  // Untreated clinical infection (U -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";

    human->house_ptr->cinc += 1;
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

  double dP = human->house_ptr->par_ptr->at("dP");
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

  double uB = human->house_ptr->par_ptr->at("uB");
  double uC = human->house_ptr->par_ptr->at("uC");
  double uD = human->house_ptr->par_ptr->at("uD");
  double dB = human->house_ptr->par_ptr->at("dB");
  double dC = human->house_ptr->par_ptr->at("dC");
  double dID = human->house_ptr->par_ptr->at("dID");
  double dM = human->house_ptr->par_ptr->at("dM");

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

  double b0 = human->house_ptr->par_ptr->at("b0");
  double b1 = human->house_ptr->par_ptr->at("b1");
  double IB0 = human->house_ptr->par_ptr->at("IB0");
  double kappaB = human->house_ptr->par_ptr->at("kappaB");
  double IB = human->IB;

  /* epsilon is my personal EIR today */
  human->epsilon = human->house_ptr->EIR * ((human->house_ptr->pi.at(human->id) * get_y(human)) / human->house_ptr->W);

  /* b is the probability of infection */
  double b = b0*(b1 + ((1.-b1)/(1. + std::pow((IB/IB0),kappaB))));

  /* lambda is my personal force of infection (my hazard) */
  human->lambda = human->epsilon * b;
};

/* phi */
void update_phi(human_ptr& human){

  double phi0 = human->house_ptr->par_ptr->at("phi0");
  double phi1 = human->house_ptr->par_ptr->at("phi1");
  double IC0 = human->house_ptr->par_ptr->at("IC0");
  double kappaC = human->house_ptr->par_ptr->at("kappaC");

  double ICA = human->ICA;
  double ICM = human->ICM;

  human->phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
void update_q(human_ptr& human){

  double fD0 = human->house_ptr->par_ptr->at("fD0");
  double aD = human->house_ptr->par_ptr->at("aD");
  double gammaD = human->house_ptr->par_ptr->at("gammaD");
  double d1 = human->house_ptr->par_ptr->at("d1");
  double ID0 = human->house_ptr->par_ptr->at("ID0");
  double kappaD = human->house_ptr->par_ptr->at("kappaD");
  double alphaA = human->house_ptr->par_ptr->at("alphaA");
  double alphaU = human->house_ptr->par_ptr->at("alphaU");

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
  human->c = human->house_ptr->par_ptr->at("cT");
};

void infectiousness_D(human_ptr& human){
  human->c = human->house_ptr->par_ptr->at("cD");
};

void infectiousness_A(human_ptr& human){
  double cU = human->house_ptr->par_ptr->at("cU");
  double cD = human->house_ptr->par_ptr->at("cD");
  double gammaI = human->house_ptr->par_ptr->at("gammaI");
  human->c = cU + (cD - cU)*std::pow(human->prDetectAMic,gammaI);
};

void infectiousness_U(human_ptr& human){
  human->c = human->house_ptr->par_ptr->at("cU");
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

    double phiI = human->house_ptr->par_ptr->at("phiI");
    double rS = human->house_ptr->par_ptr->at("rIRS");
    double sS = human->house_ptr->par_ptr->at("sIRS");

    return (1. - phiI) + (phiI * (1. - rS) * sS);
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = human->house_ptr->par_ptr->at("phiB");
    double sN = human->house_ptr->par_ptr->at("sITN");

    return (1. - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = human->house_ptr->par_ptr->at("phiI");
    double rS = human->house_ptr->par_ptr->at("rIRS");
    double sS = human->house_ptr->par_ptr->at("sIRS");

    double phiB = human->house_ptr->par_ptr->at("phiB");
    double sN = human->house_ptr->par_ptr->at("sITN");

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

    double phiI = human->house_ptr->par_ptr->at("phiI");
    double rS = human->house_ptr->par_ptr->at("rIRS");

    return (1.0 - phiI) + (phiI * (1. - rS));
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = human->house_ptr->par_ptr->at("phiB");
    double sN = human->house_ptr->par_ptr->at("sITN");

    return (1.0 - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = human->house_ptr->par_ptr->at("phiI");
    double rS = human->house_ptr->par_ptr->at("rIRS");

    double phiB = human->house_ptr->par_ptr->at("phiB");
    double sN = human->house_ptr->par_ptr->at("sITN");

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

    double phiI = human->house_ptr->par_ptr->at("phiI");
    double rS = human->house_ptr->par_ptr->at("rIRS");

    return phiI * rS;
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = human->house_ptr->par_ptr->at("phiB");
    double rN = human->house_ptr->par_ptr->at("rN");

    return phiB * rN;
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = human->house_ptr->par_ptr->at("phiI");
    double rS = human->house_ptr->par_ptr->at("rIRS");

    double phiB = human->house_ptr->par_ptr->at("phiB");
    double rN = human->house_ptr->par_ptr->at("rN");

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

  double a0 = human->house_ptr->par_ptr->at("a0");
  double rho = human->house_ptr->par_ptr->at("rho");
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

  double a0 = human->house_ptr->par_ptr->at("a0");
  double rho = human->house_ptr->par_ptr->at("rho");
  double pi = human->zeta * (1. - rho * std::exp(-human->age/a0));

  // update the hash table
  human->house_ptr->pi.at(human->id) = pi;

};


/* ################################################################################
#   Humans: daily update
################################################################################ */

void one_day_update_human(human_ptr& human, const int tnow){

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
    update_interventions_human(human,tnow);

  }

};


/* ################################################################################
#   interventions
################################################################################ */

// called before exiting daily update; check if interventions expire
void update_interventions_human(human_ptr& human, const int tnow){

  if(human->ITN && tnow >= human->ITN_decay){
    human->ITN = false;
  }

};

void give_ITN(human_ptr& human, const int tnow){

  double ITN_decay = human->house_ptr->par_ptr->at("ITN_decay");

  human->ITN = true;
  human->ITN_given = tnow;
  human->ITN_decay = tnow + R::rgeom(ITN_decay);

};


// /* ################################################################################
// #   track output
// ################################################################################ */
//
// Rcpp::NumericVector track_transmission(human const* const human_ptr){
//   return Rcpp::NumericVector::create(
//     Rcpp::Named("epsilon") = human_ptr->epsilon,
//     Rcpp::Named("lambda") = human_ptr->lambda,
//     Rcpp::Named("phi") = human_ptr->phi,
//     Rcpp::Named("prDetectAMic") = human_ptr->prDetectAMic,
//     Rcpp::Named("prDetectAPCR") = human_ptr->prDetectAPCR,
//     Rcpp::Named("prDetectUPCR") = human_ptr->prDetectUPCR,
//     Rcpp::Named("c") = human_ptr->c
//   );
// };
//
// Rcpp::NumericVector track_immunity(human const* const human_ptr){
//   return Rcpp::NumericVector::create(
//     Rcpp::Named("IB") = human_ptr->IB,
//     Rcpp::Named("ID") = human_ptr->ID,
//     Rcpp::Named("ICA") = human_ptr->ICA,
//     Rcpp::Named("ICM") = human_ptr->ICM
//   );
// };
//
// Rcpp::List human_2list(human const* const human_ptr){
//   return Rcpp::List::create(
//     Rcpp::Named("id") = human_ptr->id,
//     Rcpp::Named("age") = human_ptr->age,
//     Rcpp::Named("zeta") = human_ptr->zeta,
//     Rcpp::Named("IB") = human_ptr->IB,
//     Rcpp::Named("ID") = human_ptr->ID,
//     Rcpp::Named("ICA") = human_ptr->ICA,
//     Rcpp::Named("ICM") = human_ptr->ICM,
//     Rcpp::Named("epsilon") = human_ptr->epsilon,
//     Rcpp::Named("lambda") = human_ptr->lambda,
//     Rcpp::Named("phi") = human_ptr->phi,
//     Rcpp::Named("prDetectAMic") = human_ptr->prDetectAMic,
//     Rcpp::Named("prDetectAPCR") = human_ptr->prDetectAPCR,
//     Rcpp::Named("prDetectUPCR") = human_ptr->prDetectUPCR,
//     Rcpp::Named("c") = human_ptr->c,
//     Rcpp::Named("state") = human_ptr->state
//   );
// };
