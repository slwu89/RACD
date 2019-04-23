/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Marshall Lab
 #  February 2019
 #
 #  Mosquito Habitat Class Implementation
*/

#include "RACD-Mosquito.hpp"
#include "RACD-Village.hpp"

/* ################################################################################
 * basic methods needed
################################################################################ */

/* constructor */
mosquito_habitat::mosquito_habitat(const int EL_, const int LL_, const int PL_, const int SV_, const int EV_, const int IV_, const double K_,
  const Rcpp::NumericVector& psiR, const Rcpp::List pars_, village* const village_ptr_) :
  psi(Rcpp::as<std::vector<double> >(psiR)), village_ptr(village_ptr_),
  EL_probs{0,0,0}, EL_transitions{0,0,0},
  LL_probs{0,0,0}, LL_transitions{0,0,0},
  PL_probs{0,0,0,0}, PL_transitions{0,0,0,0},
  SV_probs{0,0,0}, SV_transitions{0,0,0},
  EV_probs{0,0,0}, EV_transitions{0,0,0},
  IV_probs{0,0}, IV_transitions{0,0},
  EL(EL_), LL(LL_), PL(PL_), SV(SV_), EV(EV_), IV(IV_), K(K_)
{
  /* initialize biological parameters */
  Rcpp::CharacterVector pars_names = pars_.names();
  for(size_t i=0; i<pars_.length(); i++){
    pars.emplace(Rcpp::as<std::string>(pars_names[i]), Rcpp::as<double>(pars_[i]));
  }

};

/* destructor */
mosquito_habitat::~mosquito_habitat(){};


/* ################################################################################
 * send out bites to houses (EIR)
################################################################################ */

void mosquito_habitat::EIR(){

  double inf_bites = 0.0; /* total num. of infectious bites */



};


/* ################################################################################
 * stochastic Euler step
################################################################################ */

void mosquito_habitat::euler_step(const double tnow, const double dt){

  /* ########################################
  # INTERVENTION-DEPENDENT PARAMETERS
  ######################################## */

  /* MOVE THESE TO PRE-SIMULATION CALCULATIONS IN PKG */
  double delta = 1.0/(pars["tau1"]+pars["tau2"]); /* Inverse of gonotrophic cycle without ITNs/IRS */
  double e_ov = pars["beta"]*(std::exp(pars["muV"]/delta)-1.0)/pars["muV"]; /* Number of eggs per oviposition per mosquito */

  /* NO INTERVENTIONS NOW */
  double ITNcov_t = 0.0;
  double IRScov_t = 0.0;

  /* zCom: Probability of a mosquito being repelled from an ITN or IRS-treated house */
  double c0 = 1.0 - ITNcov_t - IRScov_t + ITNcov_t*IRScov_t;
  double cITN = ITNcov_t - ITNcov_t*IRScov_t;
  double cIRS = IRScov_t - ITNcov_t*IRScov_t;
  double cCom = ITNcov_t*IRScov_t;
  double rCom = pars["rIRS"] + (1.0-pars["rIRS"])*pars["rITN"];
  double sCom = (1.0-pars["rIRS"])*pars["sITN"]*pars["sIRS"];
  double zCom = pars["Q0"]*cITN*pars["phiB"]*pars["rITN"] + pars["Q0"]*cIRS*pars["phiI"]*pars["rIRS"] + pars["Q0"]*cCom*(pars["phiI"]-pars["phiB"])*pars["rIRS"] + pars["Q0"]*cCom*pars["phiB"]*rCom;

  /* deltaCom: Inverse of gonotrophic cycle length with ITNs & IRS */
  double deltaCom = 1.0/(pars["tau1"]/(1-zCom) + pars["tau2"]);

  /* wCom: Probability that a surviving mosquito succeeds in feeding during a single attempt */
  double wCom = (1.0 - pars["Q0"]) + pars["Q0"]*c0 + pars["Q0"]*cITN*(1.0-pars["phiB"]+pars["phiB"]*pars["sITN"]) + pars["Q0"]*cIRS*(1.0-pars["phiI"]+pars["phiI"]*pars["sIRS"]) + pars["Q0"]*cCom*((pars["phiI"]-pars["phiB"])*pars["sIRS"] + 1.0-pars["phiI"] + pars["phiB"]*sCom);

  /* muVCom: Female mosquito death rate in presence of ITNs & IRS */
  double p10 = std::exp(-pars["muV"]*pars["tau1"]);
  double p1Com = p10*wCom/(1.0 - zCom*p10);
  double p2 = std::exp(-pars["muV"]*pars["tau2"]);
  double pCom = std::pow(p1Com*p2,deltaCom);
  double muVCom = -std::log(pCom);

  /* betaCom: Eggs laid per day by female mosquitoes in presence of ITNs & IRS */
  double betaCom = e_ov*muVCom/(std::exp(muVCom/deltaCom) - 1.0);

  /* HBR: page 6 SI: Protocol S2 of "Reducing Plasmodium falciparum Malaria Transmission in Africa: A Model-Based Evaluation of Intervention Strategies" */
  double Q = 1.0 - ((1.0 - pars["Q0"]) / wCom); /* proportion of successful bites on humans */
  double HBR = Q*deltaCom;

  /* ########################################
  # EARLY-STAGE LARVAL INSTARS (EL)
  ######################################## */

  /* inbound oviposition to EL */
  int NV = SV + EV + IV;
  EL_new = R::rpois(betaCom * NV * dt);

  /* instantaneous hazards for EL */
  double haz_EL_mort = pars["muEL"]*(1 + ((EL+LL)/K));
  double haz_EL_2LL = 1.0 / pars["durEL"];

  /* jump probabilities */
  EL_probs[0] = std::exp(-(haz_EL_mort + haz_EL_2LL)*dt);
  EL_probs[1] = (1 - EL_probs[0])*(haz_EL_mort / (haz_EL_mort + haz_EL_2LL)); /* death */
  EL_probs[2] = (1 - EL_probs[0])*(haz_EL_2LL / (haz_EL_mort + haz_EL_2LL)); /* to late-instar */

  /* sample jumps */
  rmultinom(EL, EL_probs.data(), 3, EL_transitions.data());

  /* ########################################
  # LATE-STAGE LARVAL INSTARS (LL)
  ######################################## */

  /* instantaneous hazards for LL */
  double haz_LL_mort = pars["muLL"]*(1.0 + pars["gamma"]*((EL+LL)/K));
  double haz_LL_2PL = 1.0 / pars["durLL"];

  /* jump probabilities */
  LL_probs[0] = std::exp(-(haz_LL_mort + haz_LL_2PL)*dt);
  LL_probs[1] = (1 - LL_probs[0])*(haz_LL_mort / (haz_LL_mort + haz_LL_2PL)); /* death */
  LL_probs[2] = (1 - LL_probs[0])*(haz_LL_2PL / (haz_LL_mort + haz_LL_2PL)); /* to pupae */

  /* sample jumps */
  rmultinom(LL, LL_probs.data(), 3, LL_transitions.data());

  /* ########################################
  # PUPAE (PL)
  ######################################## */

  /* instantaneous hazards for PL */
  double haz_PL_mort = pars["muPL"];
  double haz_PL_2SV_F = (1/pars["durPL"])*0.5;
  double haz_PL_2SV_M = (1/pars["durPL"])*0.5;

  /* jump probabilities */
  PL_probs[0] = std::exp(-(haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)*dt);
  PL_probs[1] = (1 - PL_probs[0])*(haz_PL_mort / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* death */
  PL_probs[2] = (1 - PL_probs[0])*(haz_PL_2SV_F / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* to susceptible female */
  PL_probs[3] = (1 - PL_probs[0])*(haz_PL_2SV_M / (haz_PL_mort + haz_PL_2SV_F + haz_PL_2SV_M)); /* to susceptible males */

  /* sample jumps */
  rmultinom(PL, PL_probs.data(), 4, PL_transitions.data());

  /* ########################################
  # SUSCEPTIBLE VECTORS (SV)
  ######################################## */

  /* instantaneous hazards for SV */
  double haz_SV_mort =  muVCom;
  double haz_SV_inf = pars["lambdaV"];

  /* jump probabilities */
  SV_probs[0] = std::exp(-(haz_SV_mort + haz_SV_inf)*dt);
  SV_probs[1] = (1 - SV_probs[0])*(haz_SV_mort / (haz_SV_mort + haz_SV_inf)); /* death */
  SV_probs[2] = (1 - SV_probs[0])*(haz_SV_inf / (haz_SV_mort + haz_SV_inf)); /* to incubating */

  /* sample jumps */
  rmultinom(SV, SV_probs.data(), 3, SV_transitions.data());

  /* ########################################
  # INCUBATING VECTORS (EV)
  ######################################## */

  /* instantaneous hazards for EV */
  double haz_EV_mort =  muVCom;
  double haz_EV_inc = 1/pars["durEV"];

  /* jump probabilities */
  EV_probs[0] = std::exp(-(haz_EV_mort + haz_EV_inc)*dt);
  EV_probs[1] = (1 - EV_probs[0])*(haz_EV_mort / (haz_EV_mort + haz_EV_inc)); /* death */
  EV_probs[2] = (1 - EV_probs[0])*(haz_EV_inc / (haz_EV_mort + haz_EV_inc)); /* to infectious */

  /* sample jumps */
  rmultinom(EV, EV_probs.data(), 3, EV_transitions.data());

  /* ########################################
  # INFECTIOUS VECTORS (IV)
  ######################################## */

  /* instantaneous hazards for IV */
  double haz_IV_mort = muVCom;

  /* jump probabilities */
  IV_probs[0] = std::exp(-haz_IV_mort*dt);
  IV_probs[1] = (1 - IV_probs[0]);

  /* sample jumps */
  rmultinom(IV, IV_probs.data(), 2, IV_transitions.data());

  /* ########################################
  # UPDATE POPULATION
  ######################################## */

  EL = EL_transitions[0] + EL_new;
  LL = LL_transitions[0] + EL_transitions[2];
  PL = PL_transitions[0] + LL_transitions[2];
  SV = SV_transitions[0] + PL_transitions[2];
  EV = EV_transitions[0] + SV_transitions[2];
  IV = IV_transitions[0] + EV_transitions[2];

};
