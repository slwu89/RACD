/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  December 2017
 #
 #  Human Class Implementation
*/

#include "RACD-Human.hpp"
#include "RACD-House.hpp"
#include "RACD-Parameters.hpp"
#include "RACD-PRNG.hpp"


/* constructor */
human::human(const int& _humanID,
      const double& _age,
      const bool& _alive,
      house* _house_ptr,
      const double& _bitingHet,
      const double& _IB,
      const double& _ID,
      const double& _ICA,
      const double& _ICM,
      const double& _epsilon,
      const double& _lambda,
      const double& _phi,
      const double& _prDetectAMic,
      const double& _prDetectAPCR,
      const double& _prDetectUPCR,
      const std::string _state,
      const int _daysLatent
    ) : humanID(_humanID), age(_age), alive(_alive), house_ptr(_house_ptr),
        bitingHet(_bitingHet), IB(_IB), ID(_ID), ICA(_ICA), ICM(_ICM),
        epsilon(_epsilon), lambda(_lambda), phi(_phi),
        prDetectAMic(_prDetectAMic), prDetectAPCR(_prDetectAPCR), prDetectUPCR(_prDetectUPCR),
        state(_state), daysLatent(_daysLatent)
{
  #ifdef DEBUG_HPP
  std::cout << "human " << humanID << " being born at " << this << std::endl;
  #endif
};

/* destructor */
human::~human(){
  #ifdef DEBUG_HPP
  std::cout << "human " << humanID << " being killed at " << this << std::endl;
  #endif
};

/* object suicide */
void human::suicide(){
  #ifdef DEBUG_HPP
  std::cout << "human " << humanID << " suiciding at " << this << std::endl;
  #endif
  for(auto it = house_ptr->get_humans().begin(); it != house_ptr->get_humans().end(); it++){
    if(it->get() == this){
      house_ptr->get_humans().erase(it);
    }
  }
}


/* Simulation Methods */

/* mortality */
void human::mortality(){

  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

  /* mu: daily death rate as a function of mean age in years */
  double mu = RACD_Parameters::instance()->get_mu();

  if(randNum <= mu){
    // REPLACE WITH REAL OUTPUT
    std::cout << "human " << humanID << " dying from natural causes" << std::endl;
    // REPLACE WITH REAL OUTPUT
    alive = false;
    suicide();
  }

};

/* S: susceptible */
void human::S_compartment(){

  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

  /* Latent infection (S -> E):
   * If the random number is less than lambda, that individual
   * develops a latent infection (E) in the next time step.
   */
  if(randNum <= lambda){
    state = "E";
    daysLatent = 1;
    // REPLACE WITH REAL OUTPUT
    std::cout << "human " << humanID << " transitioning to E" << std::endl;
    // REPLACE WITH REAL OUTPUT
  }
};


/* E: latent period */
void human::E_compartment(){
  if(daysLatent < RACD_Parameters::instance()->get_dE()){
    daysLatent++;
  } else {

    double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

    /* fT: proportion of clinical disease cases successfully treated */
    double fT = RACD_Parameters::instance()->get_fT();

    /* Treated clinical infection (E -> T):
     * If the random number is less than phi*fT, that
     * individual develops a treated clinical infection (T)
     * in the next time step.
     */
    if(randNum <= phi*fT){
      state = "T";
      daysLatent = 0;
      // REPLACE WITH REAL OUTPUT
      std::cout << "human " << humanID << " transitioning to T" << std::endl;
      // REPLACE WITH REAL OUTPUT
    }

    /* Untreated clinical infection (E -> D):
     * If the random number is greater than phi*fT and less
     * than phi, that individual develops an untreated
     * clinical infection (D) in the next time step.
     */
     if((randNum > phi*fT) & (randNum <= phi)){
       state = "D";
       daysLatent = 0;
       // REPLACE WITH REAL OUTPUT
       std::cout << "human " << humanID << " transitioning to D" << std::endl;
       // REPLACE WITH REAL OUTPUT
     }

     /* Asymptomatic infection (E -> A):
      * If the random number is greater than phi, that
      * individual develops an asymptomatic infection (A) in
      * the next time step.
      */
      if(randNum > phi){
        state = "A";
        daysLatent = 0;
        // REPLACE WITH REAL OUTPUT
        std::cout << "human " << humanID << " transitioning to A" << std::endl;
        // REPLACE WITH REAL OUTPUT
      }

  }
};

/* T: treated clinical disease */
void human::T_compartment(){

  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

  /* dT: duration of treated clinical disease (days) */
  double dT = RACD_Parameters::instance()->get_dT();

  /* Prophylactic protection (T -> P):
   * If the random number is less than 1/dT, that individual enters the
   * phase of prophylactic protection (P) in the next time step.
  */
  if(randNum <= (1/dT)){
    state = "P";
    // REPLACE WITH REAL OUTPUT
    std::cout << "human " << humanID << " transitioning to P" << std::endl;
    // REPLACE WITH REAL OUTPUT
  }

};

/* D: untreated clinical disease */
void human::D_compartment(){

  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

  /* dD: duration of untreated clinical disease (days) */
  double dD = RACD_Parameters::instance()->get_dD();

  /* Progression from diseased to asymptomatic (D -> A):
   * If the random number is less than 1/dD, that individual enters the
   * phase of asymptomatic patent infection (A) in the next time step.
  */
  if(randNum <= (1/dD)){
    state = "A";
    // REPLACE WITH REAL OUTPUT
    std::cout << "human " << humanID << " transitioning to A" << std::endl;
    // REPLACE WITH REAL OUTPUT
  }

};

/* A: asymptomatic patent (detectable by microscopy) infection */
void human::A_compartment(){

  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = RACD_Parameters::instance()->get_fT();
  /* dA: duration of patent infection (days) */
  double dA = RACD_Parameters::instance()->get_dA();

  /* Treated clinical infection (A -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
  if(randNum <= phi*fT*lambda){
    state = "T";
    // REPLACE WITH REAL OUTPUT
    std::cout << "human " << humanID << " transitioning to T" << std::endl;
    // REPLACE WITH REAL OUTPUT
  }

  /* Untreated clinical infection (A -> D):
   * If the random number is greater than phi*fT*lambda and
   * less than phi*lambda, that individual develops an
   * untreated clinical infection (D) in the next time step.
   */
   if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
     state = "D";
     // REPLACE WITH REAL OUTPUT
     std::cout << "human " << humanID << " transitioning to D" << std::endl;
     // REPLACE WITH REAL OUTPUT
   }

   /* Progression to asymptomatic sub-patent infection (A -> U):
    * If the random number is greater than phi*lambda and less
    * than (phi*lambda + 1/dA), that individual develops an asymptomatic
    * infection (A) in the next time step.
    */
    if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1/dA)))) {
      state = "U";
      // REPLACE WITH REAL OUTPUT
      std::cout << "human " << humanID << " transitioning to U" << std::endl;
      // REPLACE WITH REAL OUTPUT
    }

};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void human::U_compartment(){

  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = RACD_Parameters::instance()->get_fT();
  /* dU: duration of sub-patent infection (days) (fitted) */
  double dU = RACD_Parameters::instance()->get_dU();

  /* Treated clinical infection (U -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
   if(randNum <= phi*fT*lambda){
     state = "T";
     // REPLACE WITH REAL OUTPUT
     std::cout << "human " << humanID << " transitioning to T" << std::endl;
     // REPLACE WITH REAL OUTPUT
   }

   /* Untreated clinical infection (U -> D):
    * If the random number is greater than phi*fT*lambda and
    * less than phi*lambda, that individual develops an
    * untreated clinical infection (D) in the next time step.
    */
    if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
      state = "D";
      // REPLACE WITH REAL OUTPUT
      std::cout << "human " << humanID << " transitioning to D" << std::endl;
      // REPLACE WITH REAL OUTPUT
    }

    /* Asymptomatic infection (U -> A):
     * If the random number is greater than phi*lambda and
     * less than lambda, that individual develops a patent
     * asymptomatic infection (A) in the next time step.
     */
     if((randNum > phi*lambda) && (randNum <= lambda)){
       state = "A";
       // REPLACE WITH REAL OUTPUT
       std::cout << "human " << humanID << " transitioning to A" << std::endl;
       // REPLACE WITH REAL OUTPUT
     }

     /* Progression to asymptomatic sub-patent infection (U -> S):
      * If the random number is greater than lambda and less
      * than (lambda + 1/dU), that individual returns to the susceptible
      * state (S) in the next time step.
      */
      if((randNum > lambda) && (randNum <= (lambda + (1/dU)))){
        state = "S";
        // REPLACE WITH REAL OUTPUT
        std::cout << "human " << humanID << " transitioning to S" << std::endl;
        // REPLACE WITH REAL OUTPUT
      }

};

/* P: protection due to chemoprophylaxis treatment */
void human::P_compartment(){

    double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

    /* dP: duration of prophylactic protection following treatment (days) */
    double dP = RACD_Parameters::instance()->get_dP();

    /* Prophylactic protection (P -> S):
     * If the random number is less than 1/dP, that individual returns to
     * the susceptible state (S) in the next time step.
     */
    if(randNum <= (1/dP)){
      state = "S";
      // REPLACE WITH REAL OUTPUT
      std::cout << "human " << humanID << " transitioning to S" << std::endl;
      // REPLACE WITH REAL OUTPUT
    }

};

/* ageing */
void human::ageing(){
  age += (1/365);
};

/* immunity */
void human::update_immunity(){

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

   double uB = RACD_Parameters::instance()->get_uB();
   double uC = RACD_Parameters::instance()->get_uC();
   double uD = RACD_Parameters::instance()->get_uD();
   double dB = RACD_Parameters::instance()->get_dB();
   double dC = RACD_Parameters::instance()->get_dC();
   double dID = RACD_Parameters::instance()->get_dID();
   double dM = RACD_Parameters::instance()->get_dM();

   IB = IB + (epsilon/(epsilon*uB + 1)) - (IB)*(1/dB);
   ICA = ICA + (lambda/(lambda*uC + 1)) - (ICA)*(1/dC);
   ICM = ICM - (ICM)*(1/dM);
   ID = ID + (lambda/(lambda*uD + 1)) - (ID)*(1/dID);

};

/* lambda */
void human::update_lambda(){

  /* Lambda (the force of infection) is calculated for each individual. It
   * varies according to age and biting heterogeneity group:
   */

   double a0 = RACD_Parameters::instance()->get_a0();
   double epsilon0 = RACD_Parameters::instance()->get_epsilon0();
   double b0 = RACD_Parameters::instance()->get_b0();
   double b1 = RACD_Parameters::instance()->get_b1();
   double rho = RACD_Parameters::instance()->get_rho();
   double IB0 = RACD_Parameters::instance()->get_IB0();
   double kappaB = RACD_Parameters::instance()->get_kappaB();
   double psi = house_ptr->get_psi();

   epsilon = epsilon0 * bitingHet * (1 - rho * exp(-age/a0)) * psi;
   double b = b0*(b1 + ((1-b1)/(1 + pow((IB/IB0),kappaB))));
   lambda = epsilon * b;

};

/* phi */
void human::update_phi(){

  /* Phi (the probability of acquiring clinical disease upon infection) is
   * also calculated for each individual. It varies according to immune
   * status:
   */

   double phi0 = RACD_Parameters::instance()->get_phi0();
   double phi1 = RACD_Parameters::instance()->get_phi1();
   double IC0 = RACD_Parameters::instance()->get_IC0();
   double kappaC = RACD_Parameters::instance()->get_kappaC();

   phi = phi0 * (phi1 + ((1 - phi1)/(1 + pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
void human::update_q(){

  double fD0 = RACD_Parameters::instance()->get_fD0();
  double aD = RACD_Parameters::instance()->get_aD();
  double gammaD = RACD_Parameters::instance()->get_gammaD();
  double d1 = RACD_Parameters::instance()->get_d1();
  double ID0 = RACD_Parameters::instance()->get_ID0();
  double kappaD = RACD_Parameters::instance()->get_kappaD();
  double alphaA = RACD_Parameters::instance()->get_alphaA();
  double alphaU = RACD_Parameters::instance()->get_alphaU();

  double fD = 1 - ((1 - fD0)/(1 + pow((age/aD),gammaD)));
  double q = d1 + ((1 - d1)/(1 + (fD*pow((ID/ID0),kappaD))*fD));

  prDetectAMic = q;
  prDetectAPCR = pow(q,alphaA);
  prDetectUPCR = pow(q,alphaU);

};
