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
        ITN(false), ITNoff(2E16),
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

  /* initialize infectiousness */
  infectiousness_funs.at(state)();

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

    /* update effectiveness of interventions */
    update_intervention(tNow);

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

  double randNum = R::runif(0.0,1.0);

  /* Latent infection (S -> E) */
  if(randNum <= lambda){
    state = "E";
    daysLatent = 1;
  }
};


/* E: latent period */
void human::E_compartment(int tNow){
  if(daysLatent < house_ptr->village_ptr->param_ptr->at("dE")){
    daysLatent++;
  } else {

    double randNum = R::runif(0.0,1.0);

    /* fT: proportion of clinical disease cases successfully treated */
    double fT = house_ptr->village_ptr->param_ptr->at("fT");

    /* Treated clinical infection (E -> T):
     * If the random number is less than phi*fT, that
     * individual develops a treated clinical infection (T)
     * in the next time step.
     */
    if(randNum <= phi*fT){
      state = "T";
      daysLatent = 0;
    }

    /* Untreated clinical infection (E -> D):
     * If the random number is greater than phi*fT and less
     * than phi, that individual develops an untreated
     * clinical infection (D) in the next time step.
     */
     if((randNum > phi*fT) & (randNum <= phi)){
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
        /* simulation  */
        state = "A";
        daysLatent = 0;
      }

  }
};

/* T: treated clinical disease */
void human::T_compartment(const int tNow){

  double randNum = R::runif(0.0,1.0);

  /* dT: duration of treated clinical disease (days) */
  double dT = house_ptr->village_ptr->param_ptr->at("dT");

  /* Prophylactic protection (T -> P):
   * If the random number is less than 1/dT, that individual enters the
   * phase of prophylactic protection (P) in the next time step.
  */
  if(randNum <= (1/dT)){
    /* simulation  */
    state = "P";
  }

};

/* D: untreated clinical disease */
void human::D_compartment(const int tNow){

  double randNum = R::runif(0.0,1.0);

  /* dD: duration of untreated clinical disease (days) */
  double dD = house_ptr->village_ptr->param_ptr->at("dD");

  /* Progression from diseased to asymptomatic (D -> A):
   * If the random number is less than 1/dD, that individual enters the
   * phase of asymptomatic patent infection (A) in the next time step.
  */
  if(randNum <= (1/dD)){
    /* simulation  */
    state = "A";
  }

};

/* A: asymptomatic patent (detectable by microscopy) infection */
void human::A_compartment(const int tNow){

  double randNum = R::runif(0.0,1.0);

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = house_ptr->village_ptr->param_ptr->at("fT");
  /* dA: duration of patent infection (days) */
  double dA = house_ptr->village_ptr->param_ptr->at("dA");

  /* Treated clinical infection (A -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
  if(randNum <= phi*fT*lambda){
    /* simulation  */
    state = "T";
  }

  /* Untreated clinical infection (A -> D):
   * If the random number is greater than phi*fT*lambda and
   * less than phi*lambda, that individual develops an
   * untreated clinical infection (D) in the next time step.
   */
   if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
     /* simulation  */
     state = "D";
   }

   /* Progression to asymptomatic sub-patent infection (A -> U):
    * If the random number is greater than phi*lambda and less
    * than (phi*lambda + 1/dA), that individual develops sub-patent asymptomatic
    * infection (U) in the next time step.
    */
    if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1/dA)))) {
      /* simulation  */
      state = "U";
    }

};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void human::U_compartment(const int tNow){

  double randNum = R::runif(0.0,1.0);

  /* fT: proportion of clinical disease cases successfully treated */
  double fT = house_ptr->village_ptr->param_ptr->at("fT");
  /* dU: duration of sub-patent infection (days) (fitted) */
  double dU = house_ptr->village_ptr->param_ptr->at("dU");

  /* Treated clinical infection (U -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
   if(randNum <= phi*fT*lambda){
     /* simulation  */
     state = "T";
   }

   /* Untreated clinical infection (U -> D):
    * If the random number is greater than phi*fT*lambda and
    * less than phi*lambda, that individual develops an
    * untreated clinical infection (D) in the next time step.
    */
    if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
      /* simulation  */
      state = "D";
    }

    /* Asymptomatic infection (U -> A):
     * If the random number is greater than phi*lambda and
     * less than lambda, that individual develops a patent
     * asymptomatic infection (A) in the next time step.
     */
     if((randNum > phi*lambda) && (randNum <= lambda)){
       /* simulation  */
       state = "A";
     }

     /* Recovery to susceptible (U -> S):
      * If the random number is greater than lambda and less
      * than (lambda + 1/dU), that individual returns to the susceptible
      * state (S) in the next time step.
      */
      if((randNum > lambda) && (randNum <= (lambda + (1.0/dU)))){
        /* simulation  */
        state = "S";
      }

};

/* P: protection due to chemoprophylaxis treatment */
void human::P_compartment(const int tNow){

  double randNum = R::runif(0.0,1.0);

    /* dP: duration of prophylactic protection following treatment (days) */
    double dP = house_ptr->village_ptr->param_ptr->at("dP");

    /* Prophylactic protection (P -> S):
     * If the random number is less than 1/dP, that individual returns to
     * the susceptible state (S) in the next time step.
     */
    if(randNum <= (1/dP)){
      /* simulation  */
      state = "S";
    }

};


/* ################################################################################
#   Immunity & Ageing
################################################################################ */

/* mortality */
void human::mortality(const int tNow){

  double randNum = R::runif(0.0,1.0);

  /* mu: daily death rate as a function of mean age in years */
  double mu = house_ptr->village_ptr->param_ptr->at("mu");

  if(randNum <= mu){
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
   double uB = house_ptr->village_ptr->param_ptr->at("uB");
   double uC = house_ptr->village_ptr->param_ptr->at("uC");
   double uD = house_ptr->village_ptr->param_ptr->at("uD");
   double dB = house_ptr->village_ptr->param_ptr->at("dB");
   double dC = house_ptr->village_ptr->param_ptr->at("dC");
   double dID = house_ptr->village_ptr->param_ptr->at("dID");
   double dM = house_ptr->village_ptr->param_ptr->at("dM");

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
   double a0 = house_ptr->village_ptr->param_ptr->at("a0");
   double epsilon0 = house_ptr->village_ptr->param_ptr->at("epsilon0");
   double b0 = house_ptr->village_ptr->param_ptr->at("b0");
   double b1 = house_ptr->village_ptr->param_ptr->at("b1");
   double rho = house_ptr->village_ptr->param_ptr->at("rho");
   double IB0 = house_ptr->village_ptr->param_ptr->at("IB0");
   double kappaB = house_ptr->village_ptr->param_ptr->at("kappaB");
   double psi = 1.0;

   epsilon = epsilon0 * bitingHet * (1 - rho * std::exp(-age/a0)) * psi;
   double b = b0*(b1 + ((1-b1)/(1 + std::pow((IB/IB0),kappaB))));
   lambda = epsilon * b;

};

/* phi */
void human::update_phi(){

  /* Phi (the probability of acquiring clinical disease upon infection) is
   * also calculated for each individual. It varies according to immune
   * status:
   */

   double phi0 = house_ptr->village_ptr->param_ptr->at("phi0");
   double phi1 = house_ptr->village_ptr->param_ptr->at("phi1");
   double IC0 = house_ptr->village_ptr->param_ptr->at("IC0");
   double kappaC = house_ptr->village_ptr->param_ptr->at("kappaC");

   phi = phi0 * (phi1 + ((1 - phi1)/(1 + std::pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
void human::update_q(){

 /* q (the probability that an asymptomatic infection is detected by
  * microscopy) is also calculated for each individual, as well as the
  * probability of detection by PCR for asymptomatic infections in states
  * A (patent) and U (subpatent). This also varies according to immune status:
  */

  double fD0 = house_ptr->village_ptr->param_ptr->at("fD0");
  double aD = house_ptr->village_ptr->param_ptr->at("aD");
  double gammaD = house_ptr->village_ptr->param_ptr->at("gammaD");
  double d1 = house_ptr->village_ptr->param_ptr->at("d1");
  double ID0 = house_ptr->village_ptr->param_ptr->at("ID0");
  double kappaD = house_ptr->village_ptr->param_ptr->at("kappaD");
  double alphaA = house_ptr->village_ptr->param_ptr->at("alphaA");
  double alphaU = house_ptr->village_ptr->param_ptr->at("alphaU");

  double fD = 1 - ((1 - fD0)/(1 + std::pow((age/aD),gammaD)));
  double q = d1 + ((1 - d1)/(1 + (fD*std::pow((ID/ID0),kappaD))*fD));

  prDetectAMic = q;
  prDetectAPCR = std::pow(q,alphaA);
  prDetectUPCR = std::pow(q,alphaU);

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
  c = house_ptr->village_ptr->param_ptr->at("cT");
};

void human::infectiousness_D(){
  c = house_ptr->village_ptr->param_ptr->at("cD");
};

void human::infectiousness_A(){
  double cU = house_ptr->village_ptr->param_ptr->at("cU");
  double cD = house_ptr->village_ptr->param_ptr->at("cD");
  double gammaI = house_ptr->village_ptr->param_ptr->at("gammaI");
  c = cU + (cD - cU)*std::pow(prDetectAMic,gammaI);
};

void human::infectiousness_U(){
  c = house_ptr->village_ptr->param_ptr->at("cU");
};

void human::infectiousness_P(){
  c = 0;
};


/* ################################################################################
#   Interventions
################################################################################ */

void human::update_intervention(const int tNow){
  if(ITN && tNow > ITNoff){
    ITN = false;
  }
}

void human::apply_ITN(){
  ITN = true;
  ITNoff = R::rexp(house_ptr->village_ptr->param_ptr->at("ITNduration"));
}


/* ################################################################################
#   mosquito/biting encounter
################################################################################ */

/* my daily recieved bites */
void human::get_bitten(const int eps){

  /* no bites (yay!) */
  if(eps == 0){
    epsilon = 0;
    lambda = 0.;
  } else {

    /* calc my b P(get inf | i get infectious bite) */
    double b0 = house_ptr->village_ptr->param_ptr->at("b0");
    double b1 = house_ptr->village_ptr->param_ptr->at("b1");
    double IB0 = house_ptr->village_ptr->param_ptr->at("IB0");
    double kappaB = house_ptr->village_ptr->param_ptr->at("kappaB");
    double b = b0*(b1 + ((1-b1)/(1 + std::pow((IB/IB0),kappaB))));

    epsilon = eps;
    lambda = std::max(eps*b,1.);
    // lambda = (int)R::rbinom((double) eps, b);
  }
}

/* individual biting weight */
double human::get_pi(){
  double rho = house_ptr->village_ptr->param_ptr->at("rho");
  double a0 = house_ptr->village_ptr->param_ptr->at("a0");
  return bitingHet * (1 - rho * std::exp(-age/a0));
};

/* probability of successful biting */
double human::get_w(){
  bool IRS = house_ptr->has_IRS();
  /* none */
  if(!ITN && !IRS){
    return 1.0;
  /* IRS only */
  } else if(IRS && !ITN){

    double phiI = house_ptr->village_ptr->param_ptr->at("phiI");
    double rS = house_ptr->village_ptr->param_ptr->at("rIRS");
    double sS = house_ptr->village_ptr->param_ptr->at("sIRS");

    return (1.0 - phiI) + (phiI * (1 - rS) * sS);
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = house_ptr->village_ptr->param_ptr->at("phiB");
    double sN = house_ptr->village_ptr->param_ptr->at("sITN");

    return (1.0 - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = house_ptr->village_ptr->param_ptr->at("phiI");
    double rS = house_ptr->village_ptr->param_ptr->at("rIRS");
    double sS = house_ptr->village_ptr->param_ptr->at("sIRS");

    double phiB = house_ptr->village_ptr->param_ptr->at("phiB");
    double sN = house_ptr->village_ptr->param_ptr->at("sITN");

    return (1.0 - phiI) + (phiB * (1 - rS) * sN * sS) + ((phiI - phiB) * (1 - rS) * sS);
  } else {
    Rcpp::stop("error: invalid combination of ITN/IRS");
  }
};

/* probability of biting */
double human::get_y(){
  bool IRS = house_ptr->has_IRS();
  /* none */
  if(!ITN && !IRS){
    return 1.0;
  /* IRS only */
  } else if(IRS && !ITN){

    double phiI = house_ptr->village_ptr->param_ptr->at("phiI");
    double rS = house_ptr->village_ptr->param_ptr->at("rIRS");

    return (1.0 - phiI) + (phiI * (1 - rS));
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = house_ptr->village_ptr->param_ptr->at("phiB");
    double sN = house_ptr->village_ptr->param_ptr->at("sITN");

    return (1.0 - phiB) + (phiB * sN);
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = house_ptr->village_ptr->param_ptr->at("phiI");
    double rS = house_ptr->village_ptr->param_ptr->at("rIRS");

    double phiB = house_ptr->village_ptr->param_ptr->at("phiB");
    double sN = house_ptr->village_ptr->param_ptr->at("sITN");

    return (1.0 - phiI) + (phiB * (1 - rS) * sN) + ((phiI - phiB) * (1 - rS) );
  } else {
    Rcpp::stop("error: invalid combination of ITN/IRS");
  }
};

/* probability of repellency*/
double human::get_z(){
  bool IRS = house_ptr->has_IRS();
  /* none */
  if(!ITN && !IRS){
    return 0.0;
  /* IRS only */
  } else if(IRS && !ITN){

    double phiI = house_ptr->village_ptr->param_ptr->at("phiI");
    double rS = house_ptr->village_ptr->param_ptr->at("rIRS");

    return phiI * rS;
  /* ITN only */
  } else if(!IRS && ITN){

    double phiB = house_ptr->village_ptr->param_ptr->at("phiB");
    double rN = house_ptr->village_ptr->param_ptr->at("rN");

    return phiB * rN;
  /* IRS and ITN */
  } else if(IRS && ITN){

    double phiI = house_ptr->village_ptr->param_ptr->at("phiI");
    double rS = house_ptr->village_ptr->param_ptr->at("rIRS");

    double phiB = house_ptr->village_ptr->param_ptr->at("phiB");
    double rN = house_ptr->village_ptr->param_ptr->at("rN");

    return (phiB * (1 - rS) * rN) + (phiI * rS);
  } else {
    Rcpp::stop("error: invalid combination of ITN/IRS");
  }
};
