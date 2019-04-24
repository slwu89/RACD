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
#include "RACD-House.hpp"
#include "RACD-Human.hpp"

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
  EL(EL_), LL(LL_), PL(PL_), SV(SV_), EV(EV_), IV(IV_), K(K_),
  EIR_out(psiR.size(),0)
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
 * feeding_cycle: calculate all M->H and H->M relevant quantities (EIR,FOI,etc)
 * and account for dependence on interventions in entomological parameters
################################################################################ */

void mosquito_habitat::feeding_cycle(const double dt){

  /* constants */
  double Q0 = village_ptr->param_ptr->at("Q0");
  double tau1 = village_ptr->param_ptr->at("tau1");
  double tau2 = village_ptr->param_ptr->at("tau2");
  double muV = village_ptr->param_ptr->at("muV");
  double eggOV = village_ptr->param_ptr->at("eggOV");

  /* quantities from all humans */
  double ww = 0.0;
  double zz = 0.0;
  double cc = 0.0;
  double psi_h = 0.0;
  size_t i = 0;
  int id = 0;
  for(auto &house : village_ptr->houses){
    psi_h = psi.at(i);
    for(auto &h : house->get_humans()){
      id = h.second->get_id();
      ww += (psi_h * house->get_pi(id) * h.second->get_w());
      zz += (psi_h * house->get_pi(id) * h.second->get_z());
      cc += (psi_h * house->get_pi(id) * h.second->get_w() * h.second->get_c());
    }
    i++;
  }

  /* P(successful feed) */
  W = (1.0 - Q0) + (Q0*ww);

  /* P(repelled w/out feed) */
  Z = Q0*zz;

  /* feeding rate */
  f = 1./((tau1/(1.0 - Z)) + tau2);

  /* survival */
  double p10 = std::exp(-muV*tau1);
  p1 = p10*W/(1.0 - Z*p10);
  p2 = std::exp(-muV*tau2);
  mu = -f*std::log(p1*p1);

  /* proportion of successful bites on humans & HBR */
  Q = 1.0 - ((1.0 - Q0)/W);
  a = f*Q;

  /* distribute EIR to houses */
  size_t nhouse = village_ptr->houses.size();
  double lambda = a*IV;
  int bites = R::rpois(lambda*dt);
  rmultinom(bites,psi.data(),nhouse,EIR_out.data());

  /* calculate FOI (h->m) */
  lambdaV = a*cc;

  /* calculate egg laying rate */
  beta = eggOV*mu/(std::exp(mu/f) - 1.0);

};


/* ################################################################################
 * stochastic Euler step
################################################################################ */

void mosquito_habitat::euler_step(const double tnow, const double dt){

  double muEL = village_ptr->param_ptr->at("muEL");
  double durEL = village_ptr->param_ptr->at("durEL");

  double gamma = village_ptr->param_ptr->at("gamma");
  double muLL = village_ptr->param_ptr->at("muLL");
  double durLL = village_ptr->param_ptr->at("durLL");

  double muPL = village_ptr->param_ptr->at("muPL");
  double durPL = village_ptr->param_ptr->at("durPL");

  double durEV = village_ptr->param_ptr->at("durEV");

  /* ########################################
  # EARLY-STAGE LARVAL INSTARS (EL)
  ######################################## */

  /* inbound oviposition to EL */
  int NV = SV + EV + IV;
  EL_new = R::rpois(beta * NV * dt);

  /* instantaneous hazards for EL */
  double haz_EL_mort = muEL*(1 + ((EL+LL)/K));
  double haz_EL_2LL = 1.0 / durEL;

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
  double haz_LL_mort = muLL*(1.0 + gamma*((EL+LL)/K));
  double haz_LL_2PL = 1.0 / durLL;

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
  double haz_PL_mort = muPL;
  double haz_PL_2SV_F = (1/durPL)*0.5;
  double haz_PL_2SV_M = (1/durPL)*0.5;

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
  double haz_SV_mort =  mu;
  double haz_SV_inf = lambdaV;

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
  double haz_EV_mort =  mu;
  double haz_EV_inc = 1/durEV;

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
  double haz_IV_mort = mu;

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
