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
 #  test human functions
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

// we dont want to hurt out fingers
using uuint = unsigned int;
static size_t tnow = 0;

// output: states
// static std::vector<size_t> state_S;
// static std::vector<size_t> state_E;
// static std::vector<size_t> state_T;
// static std::vector<size_t> state_D;
// static std::vector<size_t> state_A;
// static std::vector<size_t> state_U;
// static std::vector<size_t> state_P;

// // output: pop sizes
// static std::vector<uuint> num_All;
// static std::vector<uuint> num_2_10;
// static std::vector<uuint> num_0_5;
// static std::vector<uuint> num_5_10;
// static std::vector<uuint> num_10_15;
// static std::vector<uuint> num_15Plus;
//
// // output: clinical incidence
// static std::vector<uuint> cinc_All;
// static std::vector<uuint> cinc_2_10;
// static std::vector<uuint> cinc_0_5;
// static std::vector<uuint> cinc_5_10;
// static std::vector<uuint> cinc_10_15;
// static std::vector<uuint> cinc_15Plus;


/* ################################################################################
#   Human and human population
################################################################################ */

// we're all special snowflakes here
static uuint global_hid = 0;

// a person
typedef struct human {
  // known on initialization
  uuint         id;
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
  std::string   state;
  // other dynamic members
  uuint         days_latent;
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

void mortality(human_ptr& human, const Rcpp::NumericVector& theta){
  double randNum = R::runif(0.0,1.0);
  double mu = theta["mu"];
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
void S_compartment(human_ptr& human, const Rcpp::NumericVector& theta){

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
void E_compartment(human_ptr& human, const Rcpp::NumericVector& theta){

  int dE = Rcpp::as<double>(theta["dE"]);
  double fT = Rcpp::as<double>(theta["fT"]);

  if(human->days_latent < dE){
    human->days_latent++;
  } else {
    double randNum = R::runif(0.0,1.0);

    double phi = human->phi;

    // Treated clinical infection (E -> T)
    if(randNum <= phi*fT){
      human->state = "T";
      human->days_latent = 0;
    }

    // Untreated clinical infection (E -> D)
    if((randNum > phi*fT) && (randNum <= phi)){
      human->state = "D";
      human->days_latent = 0;
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
void T_compartment(human_ptr& human, const Rcpp::NumericVector& theta){

  int dT = Rcpp::as<double>(theta["dT"]);

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
void D_compartment(human_ptr& human, const Rcpp::NumericVector& theta){

  int dD = Rcpp::as<double>(theta["dD"]);
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
void A_compartment(human_ptr& human, const Rcpp::NumericVector& theta){

  double fT = Rcpp::as<double>(theta["fT"]);
  int dA = Rcpp::as<double>(theta["dA"]);

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (A -> T)
  if(randNum <= phi*fT*lambda){
    human->state = "T";
  }
  // Untreated clinical infection (A -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";
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
void U_compartment(human_ptr& human, const Rcpp::NumericVector& theta){

  double fT = Rcpp::as<double>(theta["fT"]);
  int dU = Rcpp::as<double>(theta["dU"]);

  double randNum = R::runif(0.0,1.0);

  double phi = human->phi;
  double lambda = human->lambda;

  // Treated clinical infection (U -> T):
  if(randNum <= phi*fT*lambda){
    human->state = "T";
  }

  // Untreated clinical infection (U -> D)
  if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
    human->state = "D";
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
void P_compartment(human_ptr& human, const Rcpp::NumericVector& theta){

  int dP = Rcpp::as<double>(theta["dP"]);
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
void update_immunity(human_ptr& human, const Rcpp::NumericVector& theta){

  double uB = Rcpp::as<double>(theta["uB"]);
  double uC = Rcpp::as<double>(theta["uC"]);
  double uD = Rcpp::as<double>(theta["uD"]);
  double dB = Rcpp::as<double>(theta["dB"]);
  double dC = Rcpp::as<double>(theta["dC"]);
  double dID = Rcpp::as<double>(theta["dID"]);
  double dM = Rcpp::as<double>(theta["dM"]);

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
void update_lambda(human_ptr& human, const Rcpp::NumericVector& theta, const double psi){

  double a0 = Rcpp::as<double>(theta["a0"]);
  double epsilon0 = Rcpp::as<double>(theta["epsilon0"]);
  double b0 = Rcpp::as<double>(theta["b0"]);
  double b1 = Rcpp::as<double>(theta["b1"]);
  double rho = Rcpp::as<double>(theta["rho"]);
  double IB0 = Rcpp::as<double>(theta["IB0"]);
  double kappaB = Rcpp::as<double>(theta["kappaB"]);

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
void update_phi(human_ptr& human,const Rcpp::NumericVector& theta){

  double phi0 = Rcpp::as<double>(theta["phi0"]);
  double phi1 = Rcpp::as<double>(theta["phi1"]);
  double IC0 = Rcpp::as<double>(theta["IC0"]);
  double kappaC = Rcpp::as<double>(theta["kappaC"]);

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
void update_q(human_ptr& human,const Rcpp::NumericVector& theta){

  double fD0 = Rcpp::as<double>(theta["fD0"]);
  double aD = Rcpp::as<double>(theta["aD"]);
  double gammaD = Rcpp::as<double>(theta["gammaD"]);
  double d1 = Rcpp::as<double>(theta["d1"]);
  double ID0 = Rcpp::as<double>(theta["ID0"]);
  double kappaD = Rcpp::as<double>(theta["kappaD"]);
  double alphaA = Rcpp::as<double>(theta["alphaA"]);
  double alphaU = Rcpp::as<double>(theta["alphaU"]);

  double ID = human->ID;
  double age = human->age;

  double fD = 1. - ((1. - fD0)/(1. + std::pow((age/aD),gammaD)));
  double q = d1 + ((1. - d1)/(1. + (fD*std::pow((ID/ID0),kappaD))*fD));

  human->prDetectAMic = q;
  human->prDetectAPCR = std::pow(q,alphaA);
  human->prDetectUPCR = std::pow(q,alphaU);

};

// put the functions in a hash table
static std::unordered_map<std::string, std::function<void(human_ptr& human, const Rcpp::NumericVector&)> > state_functions  = {
  {"S",S_compartment},
  {"E",E_compartment},
  {"T",T_compartment},
  {"D",D_compartment},
  {"A",A_compartment},
  {"U",U_compartment},
  {"P",P_compartment}
};


/* ################################################################################
#   daily update
################################################################################ */

// one day for the humans
void one_day_update(const Rcpp::NumericVector& theta, const Rcpp::NumericVector& psiHouse){

  for(auto& h : human_pop){

    // mortality
    mortality(h,theta);

    // if they survived the call of the beyond
    if(h->alive){

      // state update
      std::string state(h->state);
      state_functions.at(state)(h,theta);

      // update age
      h->age += (1./365.);

      // update immunity
      update_immunity(h,theta);

      // update lambda
      double psi = psiHouse.at(h->house);
      update_lambda(h,theta,psi);

      // update phi
      update_phi(h,theta);

      // update q
      update_q(h,theta);
    }

  };

};

// a malthusian nightmare
void one_day_births(const Rcpp::NumericVector& theta, const Rcpp::NumericVector& psiHouse){

  size_t nn = human_pop.size();
  double mu = Rcpp::as<double>(theta["mu"]);

  size_t numbirth = (size_t)R::rbinom((double) nn, mu);

  if(numbirth > 0){

    double ICA18_22 = mean_ICA18_22();

    double sigma2 = Rcpp::as<double>(theta["sigma2"]);
    double epsilon0 = Rcpp::as<double>(theta["epsilon0"]);
    double rho = Rcpp::as<double>(theta["rho"]);
    double phi0 = Rcpp::as<double>(theta["phi0"]);
    double phi1 = Rcpp::as<double>(theta["phi1"]);
    double PM = Rcpp::as<double>(theta["PM"]);
    double IC0 = Rcpp::as<double>(theta["IC0"]);
    double kappaC = Rcpp::as<double>(theta["kappaC"]);
    double b0 = Rcpp::as<double>(theta["b0"]);

    for(size_t i=0; i<numbirth; i++){

      // put newborns in the smallest houses for ... reasons
      size_t smallest_house = std::distance(household_size.begin(),std::min_element(household_size.begin(),household_size.end()));
      household_size.at(smallest_house) += 1;

      // sample this person's biting heterogeneity
      double zeta = R::rlnorm(-sigma2/2., std::sqrt(sigma2));

      // their EIR
      double psi = psiHouse.at(smallest_house);
      double epsilon = epsilon0 * zeta * (1. - rho) * psi;

      // their Phi
      double phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(PM*ICA18_22/IC0,kappaC))));

      // they join the population...
      human_pop.emplace_back(std::make_unique<human>(
              0.0,
              true,
              smallest_house,
              zeta,
              0.0,
              0.0,
              0.0,
              (PM*ICA18_22),
              epsilon,
              epsilon*b0,
              phi,
              1.0,
              1.0,
              1.0,
              "S"
      ));

    }

  }

};

// bring out yer dead!
void one_day_deaths(){

  auto dead = std::remove_if(human_pop.begin(),human_pop.end(),
                             [](const human_ptr& h)
                             {return !h->alive;}
                             );
  human_pop.erase(dead,human_pop.end());

};


/* ################################################################################
#   run simulation
################################################################################ */

// hpop: list of humans
// theta: named vector of params
// psiHouse: psi values

// [[Rcpp::export]]
Rcpp::DataFrame tiny_racd_population(
  const Rcpp::List& hpop,
  const Rcpp::NumericVector& theta,
  const Rcpp::NumericVector psiHouse,
  const size_t tmax
){

  // reserve output memory
  std::vector<size_t> state_S(tmax,0);
  std::vector<size_t> state_E(tmax,0);
  std::vector<size_t> state_T(tmax,0);
  std::vector<size_t> state_D(tmax,0);
  std::vector<size_t> state_A(tmax,0);
  std::vector<size_t> state_U(tmax,0);
  std::vector<size_t> state_P(tmax,0);

  // sizes of houses
  household_size.resize(psiHouse.size(),0);

  Rcpp::Rcout << " --- initializing population in memory --- " << std::endl;
  for(size_t i=0; i<hpop.size(); i++){

    Rcpp::List human_pars = Rcpp::as<Rcpp::List>(hpop[i]);
    size_t house = Rcpp::as<size_t>(human_pars["house"]) - 1;

    household_size.at(house) += 1;

    human_pop.emplace_back(std::make_unique<human>(
      Rcpp::as<double>(human_pars["age"]),
      true,
      house,
      Rcpp::as<double>(human_pars["bitingHet"]),
      Rcpp::as<double>(human_pars["IB"]),
      Rcpp::as<double>(human_pars["ID"]),
      Rcpp::as<double>(human_pars["ICA"]),
      Rcpp::as<double>(human_pars["ICM"]),
      Rcpp::as<double>(human_pars["epsilon"]),
      Rcpp::as<double>(human_pars["lambda"]),
      Rcpp::as<double>(human_pars["phi"]),
      Rcpp::as<double>(human_pars["prDetectAMic"]),
      Rcpp::as<double>(human_pars["prDetectAPCR"]),
      Rcpp::as<double>(human_pars["prDetectUPCR"]),
      Rcpp::as<std::string>(human_pars["state"])
    ));
  }
  Rcpp::Rcout << " --- done initializing population --- " << std::endl;

  Rcpp::Rcout << " --- begin simulation --- " << std::endl;

  // main simulation loop
  Progress pb(tmax,true);
  while(tnow < tmax){
    // check for worried users
    if(tnow % 5 == 0){
      if(Progress::check_abort()){
        Rcpp::stop("user abort detected");
      }
    }

    // track output
    for(auto& h : human_pop){
      if(h->state.compare("S") == 0){
        state_S.at(tnow) += 1;
      } else if(h->state.compare("E") == 0){
        state_E.at(tnow) += 1;
      } else if(h->state.compare("T") == 0){
        state_T.at(tnow) += 1;
      } else if(h->state.compare("D") == 0){
        state_D.at(tnow) += 1;
      } else if(h->state.compare("A") == 0){
        state_A.at(tnow) += 1;
      } else if(h->state.compare("U") == 0){
        state_U.at(tnow) += 1;
      } else if(h->state.compare("P") == 0){
        state_P.at(tnow) += 1;
      } else {
        Rcpp::stop("incorrect state detected");
      }
    }

    // human simulation functions
    one_day_update(theta, psiHouse);
    one_day_births(theta, psiHouse);
    one_day_deaths();

    // increment simulation time
    pb.increment();
    tnow++;
  }

  Rcpp::Rcout << std::endl << " --- end simulation --- " << std::endl;

  // return output
  return Rcpp::DataFrame::create(
    Rcpp::Named("state_S") = Rcpp::wrap(state_S),
    Rcpp::Named("state_E") = Rcpp::wrap(state_E),
    Rcpp::Named("state_T") = Rcpp::wrap(state_T),
    Rcpp::Named("state_D") = Rcpp::wrap(state_D),
    Rcpp::Named("state_A") = Rcpp::wrap(state_A),
    Rcpp::Named("state_U") = Rcpp::wrap(state_U),
    Rcpp::Named("state_P") = Rcpp::wrap(state_P)
  );
};
