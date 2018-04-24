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
#include "RACD-Logger.hpp"


/* constructor */
human::human(const int& _humanID,
      const double& _age,
      const bool& _alive,
      const std::string& _state,
      const int& _daysLatent,
      const double& _IB,
      const double& _ID,
      const double& _ICA,
      const double& _ICM,
      const double& _bitingHet,
      const double& _epsilon,
      const double& _lambda,
      const double& _phi,
      const double& _prDetectAMic,
      const double& _prDetectAPCR,
      const double& _prDetectUPCR,
      house* _house_ptr
    ) : humanID(_humanID), age(_age), alive(_alive), state(_state), daysLatent(_daysLatent),
        IB(_IB), ID(_ID), ICA(_ICA), ICM(_ICM),
        bitingHet(_bitingHet), epsilon(_epsilon), lambda(_lambda), phi(_phi),
        prDetectAMic(_prDetectAMic), prDetectAPCR(_prDetectAPCR), prDetectUPCR(_prDetectUPCR),
        house_ptr(_house_ptr)
{
  #ifdef DEBUG_RACD
  std::cout << "human " << humanID << " being born at " << this << std::endl;
  #endif
};

/* destructor */
human::~human(){
  #ifdef DEBUG_RACD
  std::cout << "human " << humanID << " being killed at " << this << std::endl;
  #endif
};


/* Simulation Methods */

/* daily simulation */
void human::one_day(const int& tNow){

  /* daily mortality */
  mortality(tNow);

  if(alive){

    /* compartment transitions */
    if(state.compare("S")==0){
      S_compartment(tNow);
    } else if(state.compare("E")==0){
      E_compartment(tNow);
    } else if(state.compare("T")==0){
      T_compartment(tNow);
    } else if(state.compare("D")==0){
      D_compartment(tNow);
    } else if(state.compare("A")==0){
      A_compartment(tNow);
    } else if(state.compare("U")==0){
      U_compartment(tNow);
    } else if(state.compare("P")==0){
      P_compartment(tNow);
    } else {
      Rcpp::stop("unrecognized human state");
    }

    /* ageing */
    ageing();

    /* update immunity */
    update_immunity();
    update_lambda();
    update_phi();
    update_q();
  }

};

/* mortality */
void human::mortality(const int& tNow){

  double randNum = prng::instance().get_runif();

  /* mu: daily death rate as a function of mean age in years */
  double mu = RACD_Parameters::instance().get_mu();

  if(randNum <= mu){

    /* logging */
    std::string out = std::to_string(humanID) + ",Death," + std::to_string(tNow) + "," + std::to_string(age);
    logger::instance().log_trans(out);

    /* simulation */
    alive = false;
  }

};

/* S: susceptible */
void human::S_compartment(const int& tNow){

  double randNum = prng::instance().get_runif();

  /* Latent infection (S -> E):
   * If the random number is less than lambda, that individual
   * develops a latent infection (E) in the next time step.
   */
  if(randNum <= lambda){

    /* logging */
    std::string out = std::to_string(humanID) + ",E," + std::to_string(tNow) + "," + std::to_string(age);
    logger::instance().log_trans(out);

    /* simulation */
    state = "E";
    daysLatent = 1;
  }
};


/* E: latent period */
void human::E_compartment(const int& tNow){
  if(daysLatent < RACD_Parameters::instance().get_dE()){
    daysLatent++;
  } else {

    double randNum = prng::instance().get_runif();

    /* fT: proportion of clinical disease cases successfully treated */
    double fT = RACD_Parameters::instance().get_fT();

    /* Treated clinical infection (E -> T):
     * If the random number is less than phi*fT, that
     * individual develops a treated clinical infection (T)
     * in the next time step.
     */
    if(randNum <= phi*fT){

      /* logging */
      std::string out = std::to_string(humanID) + ",T," + std::to_string(tNow) + "," + std::to_string(age);
      logger::instance().log_trans(out);

      /* simulation  */
      state = "T";
      daysLatent = 0;
    }

    /* Untreated clinical infection (E -> D):
     * If the random number is greater than phi*fT and less
     * than phi, that individual develops an untreated
     * clinical infection (D) in the next time step.
     */
     if((randNum > phi*fT) & (randNum <= phi)){

       /* logging */
       std::string out = std::to_string(humanID) + ",D," + std::to_string(tNow) + "," + std::to_string(age);
       logger::instance().log_trans(out);

       /* simulation  */
       state = "D";
       daysLatent = 0;
     }

     /* Asymptomatic infection (E -> A):
      * If the random number is greater than phi, that
      * individual develops an asymptomatic infection (A) in
      * the next time step.
      */
      if(randNum > phi){

        /* logging */
        std::string out = std::to_string(humanID) + ",A," + std::to_string(tNow) + "," + std::to_string(age);
        logger::instance().log_trans(out);

        /* simulation  */
        state = "A";
        daysLatent = 0;
      }

  }
};

/* T: treated clinical disease */
void human::T_compartment(const int& tNow){

  double randNum = prng::instance().get_runif();

  /* dT: duration of treated clinical disease (days) */
  double dT = RACD_Parameters::instance().get_dT();

  /* Prophylactic protection (T -> P):
   * If the random number is less than 1/dT, that individual enters the
   * phase of prophylactic protection (P) in the next time step.
  */
  if(randNum <= (1/dT)){

    /* logging */
    std::string out = std::to_string(humanID) + ",P," + std::to_string(tNow) + "," + std::to_string(age);
    logger::instance().log_trans(out);

    /* simulation  */
    state = "P";
  }

};

/* D: untreated clinical disease */
void human::D_compartment(const int& tNow){

  double randNum = prng::instance().get_runif();

  /* dD: duration of untreated clinical disease (days) */
  double dD = RACD_Parameters::instance().get_dD();

  /* Progression from diseased to asymptomatic (D -> A):
   * If the random number is less than 1/dD, that individual enters the
   * phase of asymptomatic patent infection (A) in the next time step.
  */
  if(randNum <= (1/dD)){

    /* logging */
    std::string out = std::to_string(humanID) + ",A," + std::to_string(tNow) + "," + std::to_string(age);
    logger::instance().log_trans(out);

    /* simulation  */
    state = "A";
  }

};

/* A: asymptomatic patent (detectable by microscopy) infection */
void human::A_compartment(const int& tNow){

  double randNum = prng::instance().get_runif();

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = RACD_Parameters::instance().get_fT();
  /* dA: duration of patent infection (days) */
  double dA = RACD_Parameters::instance().get_dA();

  /* Treated clinical infection (A -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
  if(randNum <= phi*fT*lambda){

    /* logging */
    std::string out = std::to_string(humanID) + ",T," + std::to_string(tNow) + "," + std::to_string(age);
    logger::instance().log_trans(out);

    /* simulation  */
    state = "T";
  }

  /* Untreated clinical infection (A -> D):
   * If the random number is greater than phi*fT*lambda and
   * less than phi*lambda, that individual develops an
   * untreated clinical infection (D) in the next time step.
   */
   if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){

     /* logging */
     std::string out = std::to_string(humanID) + ",D," + std::to_string(tNow) + "," + std::to_string(age);
     logger::instance().log_trans(out);

     /* simulation  */
     state = "D";
   }

   /* Progression to asymptomatic sub-patent infection (A -> U):
    * If the random number is greater than phi*lambda and less
    * than (phi*lambda + 1/dA), that individual develops an asymptomatic
    * infection (A) in the next time step.
    */
    if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1/dA)))) {

      /* logging */
      std::string out = std::to_string(humanID) + ",U," + std::to_string(tNow) + "," + std::to_string(age);
      logger::instance().log_trans(out);

      /* simulation  */
      state = "U";
    }

};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void human::U_compartment(const int& tNow){

  double randNum = prng::instance().get_runif();

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = RACD_Parameters::instance().get_fT();
  /* dU: duration of sub-patent infection (days) (fitted) */
  double dU = RACD_Parameters::instance().get_dU();

  /* Treated clinical infection (U -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
   if(randNum <= phi*fT*lambda){

     /* logging */
     std::string out = std::to_string(humanID) + ",T," + std::to_string(tNow) + "," + std::to_string(age);
     logger::instance().log_trans(out);

     /* simulation  */
     state = "T";
   }

   /* Untreated clinical infection (U -> D):
    * If the random number is greater than phi*fT*lambda and
    * less than phi*lambda, that individual develops an
    * untreated clinical infection (D) in the next time step.
    */
    if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){

      /* logging */
      std::string out = std::to_string(humanID) + ",D," + std::to_string(tNow) + "," + std::to_string(age);
      logger::instance().log_trans(out);

      /* simulation  */
      state = "D";
    }

    /* Asymptomatic infection (U -> A):
     * If the random number is greater than phi*lambda and
     * less than lambda, that individual develops a patent
     * asymptomatic infection (A) in the next time step.
     */
     if((randNum > phi*lambda) && (randNum <= lambda)){

       /* logging */
       std::string out = std::to_string(humanID) + ",A," + std::to_string(tNow) + "," + std::to_string(age);
       logger::instance().log_trans(out);

       /* simulation  */
       state = "A";
     }

     /* Progression to asymptomatic sub-patent infection (U -> S):
      * If the random number is greater than lambda and less
      * than (lambda + 1/dU), that individual returns to the susceptible
      * state (S) in the next time step.
      */
      if((randNum > lambda) && (randNum <= (lambda + (1/dU)))){

        /* logging */
        std::string out = std::to_string(humanID) + ",S," + std::to_string(tNow) + "," + std::to_string(age);
        logger::instance().log_trans(out);

        /* simulation  */
        state = "S";
      }

};

/* P: protection due to chemoprophylaxis treatment */
void human::P_compartment(const int& tNow){

  double randNum = prng::instance().get_runif();

    /* dP: duration of prophylactic protection following treatment (days) */
    double dP = RACD_Parameters::instance().get_dP();

    /* Prophylactic protection (P -> S):
     * If the random number is less than 1/dP, that individual returns to
     * the susceptible state (S) in the next time step.
     */
    if(randNum <= (1/dP)){

      /* logging */
      std::string out = std::to_string(humanID) + ",S," + std::to_string(tNow) + "," + std::to_string(age);
      logger::instance().log_trans(out);

      /* simulation  */
      state = "S";
    }

};

/* ageing */
void human::ageing(){
  age = age + (1.0/365.0);
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

   double uB = RACD_Parameters::instance().get_uB();
   double uC = RACD_Parameters::instance().get_uC();
   double uD = RACD_Parameters::instance().get_uD();
   double dB = RACD_Parameters::instance().get_dB();
   double dC = RACD_Parameters::instance().get_dC();
   double dID = RACD_Parameters::instance().get_dID();
   double dM = RACD_Parameters::instance().get_dM();

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

   double a0 = RACD_Parameters::instance().get_a0();
   double epsilon0 = RACD_Parameters::instance().get_epsilon0();
   double b0 = RACD_Parameters::instance().get_b0();
   double b1 = RACD_Parameters::instance().get_b1();
   double rho = RACD_Parameters::instance().get_rho();
   double IB0 = RACD_Parameters::instance().get_IB0();
   double kappaB = RACD_Parameters::instance().get_kappaB();
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

   double phi0 = RACD_Parameters::instance().get_phi0();
   double phi1 = RACD_Parameters::instance().get_phi1();
   double IC0 = RACD_Parameters::instance().get_IC0();
   double kappaC = RACD_Parameters::instance().get_kappaC();

   phi = phi0 * (phi1 + ((1 - phi1)/(1 + pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
void human::update_q(){

 /* q (the probability that an asymptomatic infection is detected by
  * microscopy) is also calculated for each individual, as well as the
  * probability of detection by PCR for asymptomatic infections in states
  * A (patent) and U (subpatent). This also varies according to immune status:
  */

  double fD0 = RACD_Parameters::instance().get_fD0();
  double aD = RACD_Parameters::instance().get_aD();
  double gammaD = RACD_Parameters::instance().get_gammaD();
  double d1 = RACD_Parameters::instance().get_d1();
  double ID0 = RACD_Parameters::instance().get_ID0();
  double kappaD = RACD_Parameters::instance().get_kappaD();
  double alphaA = RACD_Parameters::instance().get_alphaA();
  double alphaU = RACD_Parameters::instance().get_alphaU();

  double fD = 1 - ((1 - fD0)/(1 + pow((age/aD),gammaD)));
  double q = d1 + ((1 - d1)/(1 + (fD*pow((ID/ID0),kappaD))*fD));

  prDetectAMic = q;
  prDetectAPCR = pow(q,alphaA);
  prDetectUPCR = pow(q,alphaU);

};
