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

#include "mosquito.hpp"

// [[Rcpp::plugins(cpp14)]]


/* --------------------------------------------------------------------------------
#   mosquito ctor/dtor
-------------------------------------------------------------------------------- */

/* constructor & destructor */
mosquitos::mosquitos(const int EL_, const int LL_, const int PL_, const int SV_, const int EV_, const int IV_, const double K_,
          const double lambdaV_) :
  EL_probs{0,0,0}, EL_transitions{0,0,0},
  LL_probs{0,0,0}, LL_transitions{0,0,0},
  PL_probs{0,0,0,0}, PL_transitions{0,0,0,0},
  SV_probs{0,0,0}, SV_transitions{0,0,0},
  EV_probs{0,0,0}, EV_transitions{0,0,0},
  IV_probs{0,0}, IV_transitions{0,0},
  EL(EL_), LL(LL_), PL(PL_), SV(SV_), EV(EV_), IV(IV_), K(K_),
  lambdaV(lambdaV_)
{
  Rcpp::Rcout << "mosquitos ctor called at " << this << "\n";
};

mosquitos::~mosquitos(){
  Rcpp::Rcout << "mosquitos dtor called at " << this << "\n";
};

/* --------------------------------------------------------------------------------
#   initialize mosquito and return to R
-------------------------------------------------------------------------------- */

// [[Rcpp::export]]
Rcpp::XPtr<mosquitos> init_mosquitos(
  const int EL_,
  const int LL_,
  const int PL_,
  const int SV_,
  const int EV_,
  const int IV_,
  const double K_,
  const double lambdaV_
){
  mosquitos* m_ptr = new mosquitos(EL_,LL_,PL_,SV_,EV_,IV_,K_,lambdaV_);

  // return the external pointer
  Rcpp::XPtr<mosquitos> p(m_ptr);
  return p;
};


/* --------------------------------------------------------------------------------
#   daily updates
-------------------------------------------------------------------------------- */

// feeding cycle needs to return the vector of EIR for houses;
// bite_probs: vector of WW,ZZ,CC
// psi: vector of probabilities to distribute bites to houses
// [[Rcpp::export]]
std::vector<double> feeding_cycle(SEXP mosy,
  const std::vector<double>& bite_probs,
  const std::vector<double>& psi,
  const Rcpp::NumericVector& parameters
){

  // grab the pointer to mosquitos
  Rcpp::XPtr<mosquitos> mosy_ptr(mosy);

  // fixed daily time step
  const double dt = 1.;

  /* constants */
  double Q0 = parameters["Q0"];
  double tau1 = parameters["tau1"];
  double tau2 = parameters["tau2"];
  double muV = parameters["muV"];
  double eggOV = parameters["eggOV"];

  // biting & human -> mosquito transmission
  double WW = bite_probs[0];
  double ZZ = bite_probs[1];
  double CC = bite_probs[2];
  
  /* P(successful feed) */
  mosy_ptr->W = (1.0 - Q0) + (Q0*WW);

  /* P(repelled w/out feed) */
  mosy_ptr->Z = Q0*ZZ;

  /* feeding rate */
  mosy_ptr->f = 1. /(tau1/(1. - mosy_ptr->Z) + tau2);

  /* survival */
  double p10 = std::exp(-muV*tau1);
  if(mosy_ptr->Z > (1. -  p10*mosy_ptr->W)/p10){
    std::string msg = "error: Z " + std::to_string(mosy_ptr->Z) + " is greater than 1 - p10*W / p10: " + std::to_string((1. -  p10*mosy_ptr->W)/p10) + ", p1 will have negative probabilities";
    Rcpp::stop(msg);
  }
  mosy_ptr->p1 = p10*mosy_ptr->W/(1.0 - mosy_ptr->Z*p10);
  mosy_ptr->p2 = std::exp(-muV*tau2);
  mosy_ptr->mu = -mosy_ptr->f*std::log(mosy_ptr->p1 * mosy_ptr->p2);

  /* proportion of successful bites on humans & HBR */
  mosy_ptr->Q = 1.-(1.-Q0)/mosy_ptr->W;
  mosy_ptr->a = mosy_ptr->f*mosy_ptr->Q;

  /* calculate egg laying rate */
  mosy_ptr->beta = eggOV*mosy_ptr->mu/(std::exp(mosy_ptr->mu/mosy_ptr->f) - 1.0);

  /* calculate FOI (h->m) */
  mosy_ptr->lambdaV = mosy_ptr->a*CC;

  /* calculate EIR on houses */
  double bites = mosy_ptr->a * (double)mosy_ptr->IV * dt;
  std::vector<double> EIR(psi.size());
  for(size_t i=0; i<psi.size(); i++){
    EIR[i] = psi[i] * bites;
  }
  return EIR;
};

// [[Rcpp::export]]
void euler_step(
  SEXP mosy,
  const Rcpp::NumericVector& parameters
){

  // grab the pointer to mosquitos
  Rcpp::XPtr<mosquitos> mosy_ptr(mosy);

  const double dt = 1.0;

  double muEL = parameters["muEL"];
  double durEL = parameters["durEL"];

  double gamma = parameters["gamma"];
  double muLL = parameters["muLL"];
  double durLL = parameters["durLL"];

  double muPL = parameters["muPL"];
  double durPL = parameters["durPL"];

  double durEV = parameters["durEV"];

  /* ########################################
  # EARLY-STAGE LARVAL INSTARS (EL)
  ######################################## */

  /* inbound oviposition to EL */
  int NV = mosy_ptr->SV + mosy_ptr->EV + mosy_ptr->IV;
  mosy_ptr->EL_new = R::rpois(mosy_ptr->beta * NV * dt);

  /* instantaneous hazards for EL */
  double haz_EL_mort = muEL*(1. + ((mosy_ptr->EL+mosy_ptr->LL)/mosy_ptr->K));
  double haz_EL_2LL = 1.0 / durEL;

  /* jump probabilities */
  mosy_ptr->EL_probs[0] = std::exp(-(haz_EL_mort + haz_EL_2LL)*dt);
  mosy_ptr->EL_probs[1] = (1. - mosy_ptr->EL_probs[0])*(haz_EL_mort / (haz_EL_mort + haz_EL_2LL)); /* death */
  mosy_ptr->EL_probs[2] = (1. - mosy_ptr->EL_probs[0])*(haz_EL_2LL / (haz_EL_mort + haz_EL_2LL)); /* to late-instar */

  /* sample jumps */
  rmultinom(mosy_ptr->EL, mosy_ptr->EL_probs.data(), 3, mosy_ptr->EL_transitions.data());

  /* ########################################
  # LATE-STAGE LARVAL INSTARS (LL)
  ######################################## */

  /* instantaneous hazards for LL */
  double haz_LL_mort = muLL*(1. + gamma*((mosy_ptr->EL+mosy_ptr->LL)/mosy_ptr->K));
  double haz_LL_2PL = 1.0 / durLL;

  /* jump probabilities */
  mosy_ptr->LL_probs[0] = std::exp(-(haz_LL_mort + haz_LL_2PL)*dt);
  mosy_ptr->LL_probs[1] = (1. - mosy_ptr->LL_probs[0])*(haz_LL_mort / (haz_LL_mort + haz_LL_2PL)); /* death */
  mosy_ptr->LL_probs[2] = (1. - mosy_ptr->LL_probs[0])*(haz_LL_2PL / (haz_LL_mort + haz_LL_2PL)); /* to pupae */

  /* sample jumps */
  rmultinom(mosy_ptr->LL, mosy_ptr->LL_probs.data(), 3, mosy_ptr->LL_transitions.data());

  /* ########################################
  # PUPAE (PL)
  ######################################## */

  /* instantaneous hazards for PL */
  double haz_PL_mort = muPL;
  double haz_PL_2SV_F = (1./durPL)*0.5;
  double haz_PL_2SV_M = (1./durPL)*0.5;

  /* jump probabilities */
  mosy_ptr->PL_probs[0] = std::exp(-(haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)*dt);
  mosy_ptr->PL_probs[1] = (1. - mosy_ptr->PL_probs[0])*(haz_PL_mort / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* death */
  mosy_ptr->PL_probs[2] = (1. - mosy_ptr->PL_probs[0])*(haz_PL_2SV_F / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* to susceptible female */
  mosy_ptr->PL_probs[3] = (1. - mosy_ptr->PL_probs[0])*(haz_PL_2SV_M / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* to susceptible males */

  /* sample jumps */
  rmultinom(mosy_ptr->PL, mosy_ptr->PL_probs.data(), 4, mosy_ptr->PL_transitions.data());

  /* ########################################
  # SUSCEPTIBLE VECTORS (SV)
  ######################################## */

  /* instantaneous hazards for SV */
  double haz_SV_mort =  mosy_ptr->mu;
  double haz_SV_inf = mosy_ptr->lambdaV;

  /* jump probabilities */
  mosy_ptr->SV_probs[0] = std::exp(-(haz_SV_mort + haz_SV_inf)*dt);
  mosy_ptr->SV_probs[1] = (1. - mosy_ptr->SV_probs[0])*(haz_SV_mort / (haz_SV_mort + haz_SV_inf)); /* death */
  mosy_ptr->SV_probs[2] = (1. - mosy_ptr->SV_probs[0])*(haz_SV_inf / (haz_SV_mort + haz_SV_inf)); /* to incubating */

  /* sample jumps */
  rmultinom(mosy_ptr->SV, mosy_ptr->SV_probs.data(), 3, mosy_ptr->SV_transitions.data());

  /* ########################################
  # INCUBATING VECTORS (EV)
  ######################################## */

  /* instantaneous hazards for EV */
  double haz_EV_mort =  mosy_ptr->mu;
  double haz_EV_inc = 1./durEV;

  /* jump probabilities */
  mosy_ptr->EV_probs[0] = std::exp(-(haz_EV_mort + haz_EV_inc)*dt);
  mosy_ptr->EV_probs[1] = (1. - mosy_ptr->EV_probs[0])*(haz_EV_mort / (haz_EV_mort + haz_EV_inc)); /* death */
  mosy_ptr->EV_probs[2] = (1. - mosy_ptr->EV_probs[0])*(haz_EV_inc / (haz_EV_mort + haz_EV_inc)); /* to infectious */

  /* sample jumps */
  rmultinom(mosy_ptr->EV, mosy_ptr->EV_probs.data(), 3, mosy_ptr->EV_transitions.data());

  /* ########################################
  # INFECTIOUS VECTORS (IV)
  ######################################## */

  /* instantaneous hazards for IV */
  double haz_IV_mort = mosy_ptr->mu;

  /* jump probabilities */
  mosy_ptr->IV_probs[0] = std::exp(-haz_IV_mort*dt);
  mosy_ptr->IV_probs[1] = (1. - mosy_ptr->IV_probs[0]);

  /* sample jumps */
  rmultinom(mosy_ptr->IV, mosy_ptr->IV_probs.data(), 2, mosy_ptr->IV_transitions.data());

  /* ########################################
  # UPDATE POPULATION
  ######################################## */

  mosy_ptr->EL = mosy_ptr->EL_transitions[0] + mosy_ptr->EL_new;
  mosy_ptr->LL = mosy_ptr->LL_transitions[0] + mosy_ptr->EL_transitions[2];
  mosy_ptr->PL = mosy_ptr->PL_transitions[0] + mosy_ptr->LL_transitions[2];
  mosy_ptr->SV = mosy_ptr->SV_transitions[0] + mosy_ptr->PL_transitions[2];
  mosy_ptr->EV = mosy_ptr->EV_transitions[0] + mosy_ptr->SV_transitions[2];
  mosy_ptr->IV = mosy_ptr->IV_transitions[0] + mosy_ptr->EV_transitions[2];

};


/* --------------------------------------------------------------------------------
#   track output
-------------------------------------------------------------------------------- */

// [[Rcpp::export]]
std::vector<int> track_mosquito(SEXP mosy){

  // grab the pointer to mosquitos
  Rcpp::XPtr<mosquitos> mosy_ptr(mosy);

  std::vector<int> out{mosy_ptr->SV,mosy_ptr->EV,mosy_ptr->IV};
  return out;
};

// [[Rcpp::export]]
double           track_lambdaV(SEXP mosy){

  // grab the pointer to mosquitos
  Rcpp::XPtr<mosquitos> mosy_ptr(mosy);

  return mosy_ptr->lambdaV;
};

// [[Rcpp::export]]
Rcpp::List       mosquito_2list(SEXP mosy){

  // grab the pointer to mosquitos
  Rcpp::XPtr<mosquitos> mosy_ptr(mosy);

  // probs/transitions
  Rcpp::List prob_trans = Rcpp::List::create(
    Rcpp::Named("EL_new") = mosy_ptr->EL_new,
    Rcpp::Named("EL_probs") = mosy_ptr->EL_probs,
    Rcpp::Named("EL_transitions") = mosy_ptr->EL_transitions,
    Rcpp::Named("LL_probs") = mosy_ptr->LL_probs,
    Rcpp::Named("LL_transitions") = mosy_ptr->LL_transitions,
    Rcpp::Named("PL_probs") = mosy_ptr->PL_probs,
    Rcpp::Named("PL_transitions") = mosy_ptr->PL_transitions,
    Rcpp::Named("SV_probs") = mosy_ptr->SV_probs,
    Rcpp::Named("SV_transitions") = mosy_ptr->SV_transitions,
    Rcpp::Named("EV_probs") = mosy_ptr->EV_probs,
    Rcpp::Named("EV_transitions") = mosy_ptr->EV_transitions,
    Rcpp::Named("IV_probs") = mosy_ptr->IV_probs,
    Rcpp::Named("IV_transitions") = mosy_ptr->IV_transitions
  );

  /* state space */
  Rcpp::NumericVector state_space = Rcpp::NumericVector::create(
    Rcpp::Named("EL") = mosy_ptr->EL,
    Rcpp::Named("LL") = mosy_ptr->LL,
    Rcpp::Named("PL") = mosy_ptr->PL,
    Rcpp::Named("SV") = mosy_ptr->SV,
    Rcpp::Named("EV") = mosy_ptr->EV,
    Rcpp::Named("IV") = mosy_ptr->IV
  );

  Rcpp::NumericVector pars = Rcpp::NumericVector::create(
    Rcpp::Named("K") = mosy_ptr->K,
    Rcpp::Named("W") = mosy_ptr->W,
    Rcpp::Named("Z") = mosy_ptr->Z,
    Rcpp::Named("f") = mosy_ptr->f,
    Rcpp::Named("mu") = mosy_ptr->mu,
    Rcpp::Named("p1") = mosy_ptr->p1,
    Rcpp::Named("p2") = mosy_ptr->p2,
    Rcpp::Named("Q") = mosy_ptr->Q,
    Rcpp::Named("a") = mosy_ptr->a,
    Rcpp::Named("lambdaV") = mosy_ptr->lambdaV,
    Rcpp::Named("beta") = mosy_ptr->beta
  );

  return Rcpp::List::create(
    Rcpp::Named("prob_trans") = prob_trans,
    Rcpp::Named("state_space") = state_space,
    Rcpp::Named("pars") = pars
  );
};
