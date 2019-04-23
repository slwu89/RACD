/* ################################################################################
#       ____  ___   __________
#      / __ \/   | / ____/ __ \
#     / /_/ / /| |/ /   / / / /
#    / _, _/ ___ / /___/ /_/ /
#   /_/ |_/_/  |_\____/_____/
#
#   Sean Wu & John M. Marshall
#   December 2017
#
#   Human Class Implementation
################################################################################ */

#include "RACD-Human.hpp"
#include "RACD-House.hpp"
#include "RACD-Village.hpp"
#include "RACD-Parameters.hpp"
#include "RACD-PRNG.hpp"
#include "RACD-Logger.hpp"


/* constructor */
human::human(const int humanID_,
      const double age_,
      const bool alive_,
      const std::string& state_,
      const int daysLatent_,
      const double IB_,
      const double ID_,
      const double ICA_,
      const double ICM_,
      const double bitingHet_,
      const double epsilon_,
      const double lambda_,
      const double phi_,
      const double prDetectAMic_,
      const double prDetectAPCR_,
      const double prDetectUPCR_,
      house* house_ptr_
    ) : humanID(humanID_), age(age_), alive(alive_), state(state_), daysLatent(daysLatent_),
        IB(IB_), ID(ID_), ICA(ICA_), ICM(ICM_),
        bitingHet(bitingHet_), epsilon(epsilon_), lambda(lambda_), phi(phi_),
        prDetectAMic(prDetectAMic_), prDetectAPCR(prDetectAPCR_), prDetectUPCR(prDetectUPCR_),
        house_ptr(house_ptr_)
{

  /* compartment daily update functions */
  compartment_funs.emplace("S",std::bind(&human::S_compartment,this,_1));
  compartment_funs.emplace("E",std::bind(&human::E_compartment,this,_1));
  compartment_funs.emplace("T",std::bind(&human::T_compartment,this,_1));
  compartment_funs.emplace("D",std::bind(&human::D_compartment,this,_1));
  compartment_funs.emplace("A",std::bind(&human::A_compartment,this,_1));
  compartment_funs.emplace("U",std::bind(&human::U_compartment,this,_1));
  compartment_funs.emplace("P",std::bind(&human::P_compartment,this,_1));

  /* infectiousness to mosquitoes */
  infectiousness_funs.emplace("S",std::bind(&human::infectiousness_S,this));
  infectiousness_funs.emplace("E",std::bind(&human::infectiousness_E,this));
  infectiousness_funs.emplace("T",std::bind(&human::infectiousness_T,this));
  infectiousness_funs.emplace("D",std::bind(&human::infectiousness_D,this));
  infectiousness_funs.emplace("A",std::bind(&human::infectiousness_A,this));
  infectiousness_funs.emplace("U",std::bind(&human::infectiousness_U,this));
  infectiousness_funs.emplace("P",std::bind(&human::infectiousness_P,this));

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


/* ################################################################################
#   Compartment Daily Update
################################################################################ */

/* daily simulation */
void human::one_day(const int tNow){

  /* daily mortality */
  mortality(tNow);

  if(alive){

    /* compartment transitions */
    compartment_funs.at(state)(tNow);

    /* ageing */
    ageing();

    /* update immunity */
    update_immunity();
    update_lambda();
    update_phi();
    update_q();

    /* infectiousness to mosquitoes */
    infectiousness_funs.at(state)();
  }

};


/* S: susceptible */
void human::S_compartment(const int tNow){

  double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

  /* Latent infection (S -> E):
   * If the random number is less than lambda, that individual
   * develops a latent infection (E) in the next time step.
   */
  if(randNum <= lambda){

    /* logging */
    house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",E," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

    /* simulation */
    state = "E";
    daysLatent = 1;
  }
};


/* E: latent period */
void human::E_compartment(int tNow){
  if(daysLatent < house_ptr->village_ptr->param_ptr->get_dE()){
    daysLatent++;
  } else {

    double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

    /* fT: proportion of clinical disease cases successfully treated */
    double fT = house_ptr->village_ptr->param_ptr->get_fT();

    /* Treated clinical infection (E -> T):
     * If the random number is less than phi*fT, that
     * individual develops a treated clinical infection (T)
     * in the next time step.
     */
    if(randNum <= phi*fT){

      /* logging */
      house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",T," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

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
       house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",D," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

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
        house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",A," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

        /* simulation  */
        state = "A";
        daysLatent = 0;
      }

  }
};

/* T: treated clinical disease */
void human::T_compartment(const int tNow){

  double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

  /* dT: duration of treated clinical disease (days) */
  double dT = house_ptr->village_ptr->param_ptr->get_dT();

  /* Prophylactic protection (T -> P):
   * If the random number is less than 1/dT, that individual enters the
   * phase of prophylactic protection (P) in the next time step.
  */
  if(randNum <= (1/dT)){

    /* logging */
    house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",P," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

    /* simulation  */
    state = "P";
  }

};

/* D: untreated clinical disease */
void human::D_compartment(const int tNow){

  double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

  /* dD: duration of untreated clinical disease (days) */
  double dD = house_ptr->village_ptr->param_ptr->get_dD();

  /* Progression from diseased to asymptomatic (D -> A):
   * If the random number is less than 1/dD, that individual enters the
   * phase of asymptomatic patent infection (A) in the next time step.
  */
  if(randNum <= (1/dD)){

    /* logging */
    house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",A," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

    /* simulation  */
    state = "A";
  }

};

/* A: asymptomatic patent (detectable by microscopy) infection */
void human::A_compartment(const int tNow){

  double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = house_ptr->village_ptr->param_ptr->get_fT();
  /* dA: duration of patent infection (days) */
  double dA = house_ptr->village_ptr->param_ptr->get_dA();

  /* Treated clinical infection (A -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
  if(randNum <= phi*fT*lambda){

    /* logging */
    house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",T," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

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
     house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",D," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

     /* simulation  */
     state = "D";
   }

   /* Progression to asymptomatic sub-patent infection (A -> U):
    * If the random number is greater than phi*lambda and less
    * than (phi*lambda + 1/dA), that individual develops sub-patent asymptomatic
    * infection (U) in the next time step.
    */
    if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1/dA)))) {

      /* logging */
      house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",U," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

      /* simulation  */
      state = "U";
    }

};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void human::U_compartment(const int tNow){

  double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = house_ptr->village_ptr->param_ptr->get_fT();
  /* dU: duration of sub-patent infection (days) (fitted) */
  double dU = house_ptr->village_ptr->param_ptr->get_dU();

  /* Treated clinical infection (U -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
   if(randNum <= phi*fT*lambda){

     /* logging */
     house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",T," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

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
      house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",D," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

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
       house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",A," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

       /* simulation  */
       state = "A";
     }

     /* Recovery to susceptible (U -> S):
      * If the random number is greater than lambda and less
      * than (lambda + 1/dU), that individual returns to the susceptible
      * state (S) in the next time step.
      */
      if((randNum > lambda) && (randNum <= (lambda + (1/dU)))){

        /* logging */
        house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",S," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

        /* simulation  */
        state = "S";
      }

};

/* P: protection due to chemoprophylaxis treatment */
void human::P_compartment(const int tNow){

  double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

    /* dP: duration of prophylactic protection following treatment (days) */
    double dP = house_ptr->village_ptr->param_ptr->get_dP();

    /* Prophylactic protection (P -> S):
     * If the random number is less than 1/dP, that individual returns to
     * the susceptible state (S) in the next time step.
     */
    if(randNum <= (1/dP)){

      /* logging */
      house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",S," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

      /* simulation  */
      state = "S";
    }

};


/* ################################################################################
#   Immunity & Ageing
################################################################################ */

/* mortality */
void human::mortality(const int tNow){

  double randNum = house_ptr->village_ptr->prng_ptr->get_runif();

  /* mu: daily death rate as a function of mean age in years */
  double mu = house_ptr->village_ptr->param_ptr->get_mu();

  if(randNum <= mu){

    /* logging */
    house_ptr->village_ptr->logger_ptr->get_log() << std::to_string(humanID) << ",Death," << std::to_string(tNow) << "," <<  std::to_string(age) << "\n";

    /* simulation */
    alive = false;
  }

};

/* ageing */
void human::ageing(){
  age += 1.0/365.0;
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

   double uB = house_ptr->village_ptr->param_ptr->get_uB();
   double uC = house_ptr->village_ptr->param_ptr->get_uC();
   double uD = house_ptr->village_ptr->param_ptr->get_uD();
   double dB = house_ptr->village_ptr->param_ptr->get_dB();
   double dC = house_ptr->village_ptr->param_ptr->get_dC();
   double dID = house_ptr->village_ptr->param_ptr->get_dID();
   double dM = house_ptr->village_ptr->param_ptr->get_dM();

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

   double a0 = house_ptr->village_ptr->param_ptr->get_a0();
   double epsilon0 = house_ptr->village_ptr->param_ptr->get_epsilon0();
   double b0 = house_ptr->village_ptr->param_ptr->get_b0();
   double b1 = house_ptr->village_ptr->param_ptr->get_b1();
   double rho = house_ptr->village_ptr->param_ptr->get_rho();
   double IB0 = house_ptr->village_ptr->param_ptr->get_IB0();
   double kappaB = house_ptr->village_ptr->param_ptr->get_kappaB();
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

   double phi0 = house_ptr->village_ptr->param_ptr->get_phi0();
   double phi1 = house_ptr->village_ptr->param_ptr->get_phi1();
   double IC0 = house_ptr->village_ptr->param_ptr->get_IC0();
   double kappaC = house_ptr->village_ptr->param_ptr->get_kappaC();

   phi = phi0 * (phi1 + ((1 - phi1)/(1 + pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
void human::update_q(){

 /* q (the probability that an asymptomatic infection is detected by
  * microscopy) is also calculated for each individual, as well as the
  * probability of detection by PCR for asymptomatic infections in states
  * A (patent) and U (subpatent). This also varies according to immune status:
  */

  double fD0 = house_ptr->village_ptr->param_ptr->get_fD0();
  double aD = house_ptr->village_ptr->param_ptr->get_aD();
  double gammaD = house_ptr->village_ptr->param_ptr->get_gammaD();
  double d1 = house_ptr->village_ptr->param_ptr->get_d1();
  double ID0 = house_ptr->village_ptr->param_ptr->get_ID0();
  double kappaD = house_ptr->village_ptr->param_ptr->get_kappaD();
  double alphaA = house_ptr->village_ptr->param_ptr->get_alphaA();
  double alphaU = house_ptr->village_ptr->param_ptr->get_alphaU();

  double fD = 1 - ((1 - fD0)/(1 + pow((age/aD),gammaD)));
  double q = d1 + ((1 - d1)/(1 + (fD*pow((ID/ID0),kappaD))*fD));

  prDetectAMic = q;
  prDetectAPCR = pow(q,alphaA);
  prDetectUPCR = pow(q,alphaU);

};


/* ################################################################################
#   Infectiousness to Mosquitoes
################################################################################ */

void human::infectiousness_S(){
  c = 0;
};

void human::infectiousness_E(){
  c = 0;
};

void human::infectiousness_T(){
  c = house_ptr->village_ptr->param_ptr->get_cT();
};

void human::infectiousness_D(){
  c = house_ptr->village_ptr->param_ptr->get_cD();
};

void human::infectiousness_A(){
  double cU = house_ptr->village_ptr->param_ptr->get_cU();
  double cD = house_ptr->village_ptr->param_ptr->get_cD();
  double gammaI = house_ptr->village_ptr->param_ptr->get_gammaI();
  c = cU + (cD - cU)*pow(prDetectAMic,gammaI);
};

void human::infectiousness_U(){
  c = house_ptr->village_ptr->param_ptr->get_cU();
};

void human::infectiousness_P(){
  c = 0;
};
