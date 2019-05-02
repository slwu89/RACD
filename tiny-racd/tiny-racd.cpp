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
 #  the whole enchilada
*/

/* ################################################################################
#   Includes & auxiliary
################################################################################ */

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppProgress)]]

// for Rcpp and R's RNG because we're lazy
#include <Rcpp.h>
#include <Rmath.h>
#include <progress.hpp>

// the best debugger
#include <iostream>

// for holding the people
#include <list>

// name one C++ project that doesn't have vectors in it...go!
#include <vector>

// for smart pointers
#include <memory>

// for holding the functions
#include <unordered_map>
#include <functional>

// for strings of course
#include <string>

// to pretend we know CS
#include <algorithm>

// for da maths
#include <math.h>


/* ################################################################################
#   output and parameters/global stuff
################################################################################ */

// current simulation time
static size_t tnow = 0;

// parameters
static std::unordered_map<std::string,double> parameters;


/* ################################################################################
#   House
################################################################################ */

// a house
// all expectations are averaging out heterogeneity in the people at this house
typedef struct house {

  // data members
  size_t            id;
  double            psi;  // P(a bite going anywhere goes here)
  double            w;    // E[P(feed and survive)]
  double            y;    // E[P(feed)]
  double            z;    // E[P(repelled without feeding)]
  double            c;    // E[P(a feed here will result in a mosquito infection)]
  size_t            n;    // number of people here
  double            EIR;  // the number of bites this house gets today
  bool              IRS;  // does my house have IRS

  house(const size_t id_, const double psi_, const double w_, const double y_, const double z_, const double c_,
        const size_t n_
  );
  ~house();
} house;

// constructor
house::house(const size_t id_, const double psi_, const double w_, const double y_, const double z_, const double c_,
             const size_t n_
) :
  id(id_), psi(psi_), w(w_), y(y_), z(z_), c(c_), n(n_), EIR(0), IRS(false)
{};

// destructor
house::~house(){};

// we're men of taste, who use smart pointers
using house_ptr = std::unique_ptr<house>;

// the houses!
static std::vector<house_ptr> houses;


/* ################################################################################
#   Human and human population
################################################################################ */

// we're all special snowflakes here
static size_t global_hid = 0;

// a person
typedef struct human {
  // known on initialization
  size_t        id;
  double        age;
  bool          alive;
  size_t        house;
  double        zeta;
  double        IB;     // Pre-erythrocytic immunity (IB, reduces the probability of infection following an infectious challenge)
  double        ID;
  double        ICA;
  double        ICM;
  double        epsilon;
  double        lambda;
  double        phi;
  double        prDetectAMic;
  double        prDetectAPCR;
  double        prDetectUPCR;
  double        c;
  std::string   state;
  // other dynamic members
  size_t         days_latent;
  // constructor/destructor
  human(const double age_,
        const bool alive_,
        const size_t house_,
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
        const std::string state_);
  ~human();
} human;

// constructor
human::human(const double age_,
      const bool alive_,
      const size_t house_,
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
      const std::string state_) :
      id(global_hid),
      age(age_),
      alive(alive_),
      house(house_),
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
      c(0.0),
      state(state_),
      days_latent(0)
      {
        // after we take our id, increment for the next person!
        global_hid++;
      };

// destructor
human::~human(){};

// the smart pointer for a human
using human_ptr = std::unique_ptr<human>;

// the humans!
static std::list<human_ptr> human_pop;

// we need to know how big each household is
static std::vector<size_t> household_size(1,0);

/* shamelessly "referenced" from Knuth, The Art of Computer Programming Vol 2, section 4.2.2 */
// the mean immunity among 18-22 year olds
static inline double mean_ICA18_22(){
  double avg = 0.;
  int t = 1;
  for (auto& h : human_pop) {
    if((h->age >= 18.) && (h->age < 22.)){
      avg += (h->ICA - avg) / t;
      ++t;
    }
  }
  // weird hack if there aren't any 18-22 yr olds
  // just take avg ICA of entire pop
  if(avg < 0.001){
    avg = 0.;
    t = 1;
    for (auto& h : human_pop) {
      avg += (h->ICA - avg) / t;
      ++t;
    }
  }
  return avg;
}


/* ################################################################################
#   State transitions for our little Markov humans
################################################################################ */

void mortality(human_ptr& human){
  double randNum = R::runif(0.0,1.0);
  double mu = parameters.at("mu");
  if(randNum <= mu){
    human->alive = false;
    household_size.at(human->house) -= 1;
  }
}


/* S: susceptible */
/* Latent infection (S -> E):
 * If the random number is less than lambda, that individual
 * develops a latent infection (E) in the next time step.
 */
void S_compartment(human_ptr& human){

  double randNum = R::runif(0.0,1.0);

  if(randNum <= human->lambda){
    human->state = "E";
    human->days_latent = 0;
  }
};


/* E: latent period */
/* Treated clinical infection (E -> T):
 * If the random number is less than phi*fT, that
 * individual develops a treated clinical infection (T)
 * in the next time step.
 */
/* Untreated clinical infection (E -> D):
* If the random number is greater than phi*fT and less
* than phi, that individual develops an untreated
* clinical infection (D) in the next time step.
*/
/* Asymptomatic infection (E -> A):
 * If the random number is greater than phi, that
 * individual develops an asymptomatic infection (A) in
 * the next time step.
 */
void E_compartment(human_ptr& human){

  double dE = parameters.at("dE");
  double fT = parameters.at("fT");

  if(human->days_latent < dE){
    human->days_latent++;
  } else {
    double randNum = R::runif(0.0,1.0);

    double phi = human->phi;

    // Treated clinical infection (E -> T)
    if(randNum <= phi*fT){
      human->state = "T";
      human->days_latent = 0;

      // track_cinc(human);
    }

    // Untreated clinical infection (E -> D)
    if((randNum > phi*fT) && (randNum <= phi)){
      human->state = "D";
      human->days_latent = 0;

      // track_cinc(human);
    }

    // Asymptomatic infection (E -> A)
    if(randNum > phi){
      human->state = "A";
      human->days_latent = 0;
    }
  }
};

/* T: treated clinical disease */
/* Prophylactic protection (T -> P):
 * If the random number is less than 1/dT, that individual enters the
 * phase of prophylactic protection (P) in the next time step.
*/
void T_compartment(human_ptr& human){

  double dT = parameters.at("dT");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dT)){
    human->state = "P";
  }

};

/* D: untreated clinical disease */
/* Progression from diseased to asymptomatic (D -> A):
 * If the random number is less than 1/dD, that individual enters the
 * phase of asymptomatic patent infection (A) in the next time step.
*/
void D_compartment(human_ptr& human){

  double dD = parameters.at("dD");
  double randNum = R::runif(0.0,1.0);

  if(randNum <= (1.0/dD)){
    human->state = "A";
  }
};

/* A: asymptomatic patent (detectable by microscopy) infection */
/* Treated clinical infection (A -> T):
 * If the random number is less than phi*fT*lambda, that
 * individual develops a treated clinical infection (T) in
 * the next time step.
 */
/* Untreated clinical infection (A -> D):
* If the random number is greater than phi*fT*lambda and
* less than phi*lambda, that individual develops an
* untreated clinical infection (D) in the next time step.
*/
/* Progression to asymptomatic sub-patent infection (A -> U):
 * If the random number is greater than phi*lambda and less
 * than (phi*lambda + 1/dA), that individual develops sub-patent asymptomatic
 * infection (U) in the next time step.
 */
void A_compartment(human_ptr& human){

  double fT = parameters.at("fT");
  double dA = parameters.at("dA");

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (A -> T)
  if(randNum <= phi*fT*lambda){
    human->state = "T";

    // track_cinc(human);
  }
  // Untreated clinical infection (A -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";

    // track_cinc(human);
  }
  // Progression to asymptomatic sub-patent infection (A -> U):
  if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1.0/dA)))) {
    human->state = "U";
  }
};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
/* Treated clinical infection (U -> T):
 * If the random number is less than phi*fT*lambda, that
 * individual develops a treated clinical infection (T) in
 * the next time step.
 */
/* Untreated clinical infection (U -> D):
* If the random number is greater than phi*fT*lambda and
* less than phi*lambda, that individual develops an
* untreated clinical infection (D) in the next time step.
*/
/* Asymptomatic infection (U -> A):
* If the random number is greater than phi*lambda and
* less than lambda, that individual develops a patent
* asymptomatic infection (A) in the next time step.
*/
/* Recovery to susceptible (U -> S):
* If the random number is greater than lambda and less
* than (lambda + 1/dU), that individual returns to the susceptible
* state (S) in the next time step.
*/
void U_compartment(human_ptr& human){

  double fT = parameters.at("fT");
  double dU = parameters.at("dU");

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (U -> T):
  if(randNum <= phi*fT*lambda){
    human->state = "T";

    // track_cinc(human);
  }

  // Untreated clinical infection (U -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";

    // track_cinc(human);
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
/* Prophylactic protection (P -> S):
 * If the random number is less than 1/dP, that individual returns to
 * the susceptible state (S) in the next time step.
 */
void P_compartment(human_ptr& human){

  double dP = parameters.at("dP");
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

  double uB = parameters.at("uB");
  double uC = parameters.at("uC");
  double uD = parameters.at("uD");
  double dB = parameters.at("dB");
  double dC = parameters.at("dC");
  double dID = parameters.at("dID");
  double dM = parameters.at("dM");

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
/* Lambda (the force of infection) is calculated for each individual. It
 * varies according to age and biting heterogeneity group:
 */
// psi is the psi of my house
void update_lambda(human_ptr& human, const double psi){

  double a0 = parameters.at("a0");
  double epsilon0 = parameters.at("epsilon0");
  double b0 = parameters.at("b0");
  double b1 = parameters.at("b1");
  double rho = parameters.at("rho");
  double IB0 = parameters.at("IB0");
  double kappaB = parameters.at("kappaB");

  double zeta = human->zeta;
  double age = human->age;
  double IB = human->IB;

  double epsilon = epsilon0 * zeta * (1. - rho * std::exp(-age/a0)) * psi;

  human->epsilon = epsilon;
  double b = b0*(b1 + ((1.-b1)/(1. + std::pow((IB/IB0),kappaB))));
  human->lambda = epsilon * b;

};

/* phi */
/* Phi (the probability of acquiring clinical disease upon infection) is
 * also calculated for each individual. It varies according to immune
 * status:
 */
void update_phi(human_ptr& human){

  double phi0 = parameters.at("phi0");
  double phi1 = parameters.at("phi1");
  double IC0 = parameters.at("IC0");
  double kappaC = parameters.at("kappaC");

  double ICA = human->ICA;
  double ICM = human->ICM;

  human->phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */
/* q (the probability that an asymptomatic infection is detected by
 * microscopy) is also calculated for each individual, as well as the
 * probability of detection by PCR for asymptomatic infections in states
 * A (patent) and U (subpatent). This also varies according to immune status:
 */
void update_q(human_ptr& human){

  double fD0 = parameters.at("fD0");
  double aD = parameters.at("aD");
  double gammaD = parameters.at("gammaD");
  double d1 = parameters.at("d1");
  double ID0 = parameters.at("ID0");
  double kappaD = parameters.at("kappaD");
  double alphaA = parameters.at("alphaA");
  double alphaU = parameters.at("alphaU");

  double ID = human->ID;
  double age = human->age;

  double fD = 1. - ((1. - fD0)/(1. + std::pow((age/aD),gammaD)));
  double q = d1 + ((1. - d1)/(1. + (fD*std::pow((ID/ID0),kappaD))*fD));

  human->prDetectAMic = q;
  human->prDetectAPCR = std::pow(q,alphaA);
  human->prDetectUPCR = std::pow(q,alphaU);

};

// put the functions in a hash table
static std::unordered_map<std::string, std::function<void(human_ptr& human)> > state_functions  = {
  {"S",S_compartment},
  {"E",E_compartment},
  {"T",T_compartment},
  {"D",D_compartment},
  {"A",A_compartment},
  {"U",U_compartment},
  {"P",P_compartment}
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
  human->c = parameters.at("cT");
};

void infectiousness_D(human_ptr& human){
  human->c = parameters.at("cD");
};

void infectiousness_A(human_ptr& human){
  double cU = parameters.at("cU");
  double cD = parameters.at("cD");
  double gammaI = parameters.at("gammaI");
  human->c = cU + (cD - cU)*std::pow(human->prDetectAMic,gammaI);
};

void infectiousness_U(human_ptr& human){
  human->c = parameters.at("cU");
};

void infectiousness_P(human_ptr& human){
  human->c = 0.;
};

// put the functions in a hash table
static std::unordered_map<std::string, std::function<void(human_ptr& human)> > infectivity_functions  = {
  {"S",infectiousness_S},
  {"E",infectiousness_E},
  {"T",infectiousness_T},
  {"D",infectiousness_D},
  {"A",infectiousness_A},
  {"U",infectiousness_U},
  {"P",infectiousness_P}
};
