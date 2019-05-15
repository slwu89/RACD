/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  May 2019
 #
 #  the mosquitos
*/

// header file
#include "mosquito.hpp"
#include "globals.hpp"


/* constructor & destructor */
mosquitos::mosquitos(const int EL_, const int LL_, const int PL_, const int SV_, const int EV_, const int IV_, const double K_) :
  EL_probs{0,0,0}, EL_transitions{0,0,0},
  LL_probs{0,0,0}, LL_transitions{0,0,0},
  PL_probs{0,0,0,0}, PL_transitions{0,0,0,0},
  SV_probs{0,0,0}, SV_transitions{0,0,0},
  EV_probs{0,0,0}, EV_transitions{0,0,0},
  IV_probs{0,0}, IV_transitions{0,0},
  EL(EL_), LL(LL_), PL(PL_), SV(SV_), EV(EV_), IV(IV_), K(K_)
{};

mosquitos::~mosquitos(){};


/* ################################################################################
#   track output
################################################################################ */

void track_mosquito(const mosquito_ptr& mosy){
  mosy_S.at(tnow) = mosy->SV;
  mosy_E.at(tnow) = mosy->EV;
  mosy_I.at(tnow) = mosy->IV;
};


/* ################################################################################
#   daily updates
################################################################################ */

// calc. intervention dependent parameters and send out bites, recieve the FOI
void feeding_cycle(mosquito_ptr& mosy){

  const double dt = 1.;

  /* constants */
  double Q0 = parameters.at("Q0");
  double tau1 = parameters.at("tau1");
  double tau2 = parameters.at("tau2");
  double muV = parameters.at("muV");
  double eggOV = parameters.at("eggOV");

  /* P(successful feed) */
  mosy->W = (1.0 - Q0) + (Q0*WW);

  /* P(repelled w/out feed) */
  mosy->Z = Q0*ZZ;

  /* feeding rate */
  mosy->f = 1.0 / (tau1/(1.0 - mosy->Z) + tau2);

  /* survival */
  double p10 = std::exp(-muV*tau1);
  mosy->p1 = p10*mosy->W/(1.0 - mosy->Z*p10);
  mosy->p2 = std::exp(-muV*tau2);
  double p12 = std::pow(mosy->p1 * mosy->p2, mosy->f);
  mosy->mu = -std::log(p12); // I FUCKED UP THIS EQUATION

  /* proportion of successful bites on humans & HBR */
  mosy->Q = 1.0 - ((1.0 - Q0)/mosy->W);
  mosy->a = mosy->f*mosy->Q;

  /* calculate egg laying rate */
  mosy->beta = eggOV*mosy->mu/(std::exp(mosy->mu/mosy->f) - 1.0);

  /* calculate FOI (h->m) */
  mosy->lambdaV = mosy->a*CC;

  /* calculate EIR (m->h) */
  double bites = mosy->a * (double)mosy->IV * dt;

  // for(size_t h=0; h<GLOBAL_PSI.size(); h++){
  for(size_t h=0; h<NHOUSE; h++){
    EIR.at(h) = GLOBAL_PSI.at(h) * bites;
  }

};

// one-day update function
void euler_step(mosquito_ptr& mosy){

  const double dt = 1.0;

  double muEL = parameters.at("muEL");
  double durEL = parameters.at("durEL");

  double gamma = parameters.at("gamma");
  double muLL = parameters.at("muLL");
  double durLL = parameters.at("durLL");

  double muPL = parameters.at("muPL");
  double durPL = parameters.at("durPL");

  double durEV = parameters.at("durEV");

  /* ########################################
  # EARLY-STAGE LARVAL INSTARS (EL)
  ######################################## */

  /* inbound oviposition to EL */
  int NV = mosy->SV + mosy->EV + mosy->IV;
  mosy->EL_new = R::rpois(mosy->beta * NV * dt);

  /* instantaneous hazards for EL */
  double haz_EL_mort = muEL*(1. + ((mosy->EL+mosy->LL)/mosy->K));
  double haz_EL_2LL = 1.0 / durEL;

  /* jump probabilities */
  mosy->EL_probs[0] = std::exp(-(haz_EL_mort + haz_EL_2LL)*dt);
  mosy->EL_probs[1] = (1. - mosy->EL_probs[0])*(haz_EL_mort / (haz_EL_mort + haz_EL_2LL)); /* death */
  mosy->EL_probs[2] = (1. - mosy->EL_probs[0])*(haz_EL_2LL / (haz_EL_mort + haz_EL_2LL)); /* to late-instar */

  /* sample jumps */
  rmultinom(mosy->EL, mosy->EL_probs.data(), 3, mosy->EL_transitions.data());

  /* ########################################
  # LATE-STAGE LARVAL INSTARS (LL)
  ######################################## */

  /* instantaneous hazards for LL */
  double haz_LL_mort = muLL*(1. + gamma*((mosy->EL+mosy->LL)/mosy->K));
  double haz_LL_2PL = 1.0 / durLL;

  /* jump probabilities */
  mosy->LL_probs[0] = std::exp(-(haz_LL_mort + haz_LL_2PL)*dt);
  mosy->LL_probs[1] = (1. - mosy->LL_probs[0])*(haz_LL_mort / (haz_LL_mort + haz_LL_2PL)); /* death */
  mosy->LL_probs[2] = (1. - mosy->LL_probs[0])*(haz_LL_2PL / (haz_LL_mort + haz_LL_2PL)); /* to pupae */

  /* sample jumps */
  rmultinom(mosy->LL, mosy->LL_probs.data(), 3, mosy->LL_transitions.data());

  /* ########################################
  # PUPAE (PL)
  ######################################## */

  /* instantaneous hazards for PL */
  double haz_PL_mort = muPL;
  double haz_PL_2SV_F = (1./durPL)*0.5;
  double haz_PL_2SV_M = (1./durPL)*0.5;

  /* jump probabilities */
  mosy->PL_probs[0] = std::exp(-(haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)*dt);
  mosy->PL_probs[1] = (1. - mosy->PL_probs[0])*(haz_PL_mort / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* death */
  mosy->PL_probs[2] = (1. - mosy->PL_probs[0])*(haz_PL_2SV_F / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* to susceptible female */
  mosy->PL_probs[3] = (1. - mosy->PL_probs[0])*(haz_PL_2SV_M / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* to susceptible males */

  /* sample jumps */
  rmultinom(mosy->PL, mosy->PL_probs.data(), 4, mosy->PL_transitions.data());

  /* ########################################
  # SUSCEPTIBLE VECTORS (SV)
  ######################################## */

  /* instantaneous hazards for SV */
  double haz_SV_mort =  mosy->mu;
  double haz_SV_inf = mosy->lambdaV;

  /* jump probabilities */
  mosy->SV_probs[0] = std::exp(-(haz_SV_mort + haz_SV_inf)*dt);
  mosy->SV_probs[1] = (1. - mosy->SV_probs[0])*(haz_SV_mort / (haz_SV_mort + haz_SV_inf)); /* death */
  mosy->SV_probs[2] = (1. - mosy->SV_probs[0])*(haz_SV_inf / (haz_SV_mort + haz_SV_inf)); /* to incubating */

  /* sample jumps */
  rmultinom(mosy->SV, mosy->SV_probs.data(), 3, mosy->SV_transitions.data());

  /* ########################################
  # INCUBATING VECTORS (EV)
  ######################################## */

  /* instantaneous hazards for EV */
  double haz_EV_mort =  mosy->mu;
  double haz_EV_inc = 1./durEV;

  /* jump probabilities */
  mosy->EV_probs[0] = std::exp(-(haz_EV_mort + haz_EV_inc)*dt);
  mosy->EV_probs[1] = (1. - mosy->EV_probs[0])*(haz_EV_mort / (haz_EV_mort + haz_EV_inc)); /* death */
  mosy->EV_probs[2] = (1. - mosy->EV_probs[0])*(haz_EV_inc / (haz_EV_mort + haz_EV_inc)); /* to infectious */

  /* sample jumps */
  rmultinom(mosy->EV, mosy->EV_probs.data(), 3, mosy->EV_transitions.data());

  /* ########################################
  # INFECTIOUS VECTORS (IV)
  ######################################## */

  /* instantaneous hazards for IV */
  double haz_IV_mort = mosy->mu;

  /* jump probabilities */
  mosy->IV_probs[0] = std::exp(-haz_IV_mort*dt);
  mosy->IV_probs[1] = (1. - mosy->IV_probs[0]);

  /* sample jumps */
  rmultinom(mosy->IV, mosy->IV_probs.data(), 2, mosy->IV_transitions.data());

  /* ########################################
  # UPDATE POPULATION
  ######################################## */

  mosy->EL = mosy->EL_transitions[0] + mosy->EL_new;
  mosy->LL = mosy->LL_transitions[0] + mosy->EL_transitions[2];
  mosy->PL = mosy->PL_transitions[0] + mosy->LL_transitions[2];
  mosy->SV = mosy->SV_transitions[0] + mosy->PL_transitions[2];
  mosy->EV = mosy->EV_transitions[0] + mosy->SV_transitions[2];
  mosy->IV = mosy->IV_transitions[0] + mosy->EV_transitions[2];

};
