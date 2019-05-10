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
mosquitos::mosquitos(const size_t EL_, const size_t LL_, const size_t PL_, const size_t SV_, const size_t EV_, const size_t IV_, const double K_) :
  EL(EL_), LL(LL_), PL(PL_), SV(SV_), EV(EV_), IV(IV_), K(K_)
{};

mosquitos::~mosquitos(){};


/* ################################################################################
#   daily updates
################################################################################ */

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
  mosy->mu = -mosy->f*std::log(mosy->p1*mosy->p1);

  /* proportion of successful bites on humans & HBR */
  mosy->Q = 1.0 - ((1.0 - Q0)/mosy->W);
  mosy->a = mosy->f*mosy->Q;

  /* calculate egg laying rate */
  mosy->beta = eggOV*mosy->mu/(std::exp(mosy->mu/mosy->f) - 1.0);

  /* calculate FOI (h->m) */
  mosy->lambdaV = mosy->a*CC;

  /* calculate EIR (m->h) */
  double bites = mosy->a * (double)mosy->IV;
  for(size_t h=0; h<psi.size(); h++){
    EIR.at(h) = psi.at(h) * bites;
  }

};

void euler_step(mosquito_ptr& mosy){

};
