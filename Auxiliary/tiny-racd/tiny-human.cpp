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

// [[Rcpp::plugins(cpp14)]]

#include <Rcpp.h>
#include <Rmath.h>

#include <iostream>
#include <unordered_map>
#include <functional>

#include <math.h>


/* ################################################################################
#   State Transitions
################################################################################ */

// [[Rcpp::export]]
void mortality(Rcpp::List& human, const Rcpp::NumericVector& theta){
  double randNum = R::runif(0.0,1.0);
  double mu = theta["mu"];
  if(randNum <= mu){
    human["alive"] = false;
  }
}


/* S: susceptible */

// [[Rcpp::export]]
void S_compartment(Rcpp::List& human, const Rcpp::NumericVector& theta){

  double randNum = R::runif(0.0,1.0);
  double lambda = Rcpp::as<double>(human["lambda"]);

  /* Latent infection (S -> E):
   * If the random number is less than lambda, that individual
   * develops a latent infection (E) in the next time step.
   */
  if(randNum <= lambda){
    human["state"] = "E";
    human["daysLatent"] = 0;
  }
};


/* E: latent period */

// [[Rcpp::export]]
void E_compartment(Rcpp::List& human, const Rcpp::NumericVector& theta){

  int dE = Rcpp::as<double>(theta["dE"]);
  double fT = Rcpp::as<double>(theta["fT"]);

  int daysLatent = Rcpp::as<int>(human["daysLatent"]);
  if(daysLatent < dE){
    // Rcpp::Rcout << " daysLatent " << daysLatent << "\n";
    human["daysLatent"] = daysLatent + 1;
  } else {
    double randNum = R::runif(0.0,1.0);
    double phi = Rcpp::as<double>(human["phi"]);
    /* Treated clinical infection (E -> T):
     * If the random number is less than phi*fT, that
     * individual develops a treated clinical infection (T)
     * in the next time step.
     */
    if(randNum <= phi*fT){
      human["state"] = "T";
      human["daysLatent"] = 0;
    }
    /* Untreated clinical infection (E -> D):
     * If the random number is greater than phi*fT and less
     * than phi, that individual develops an untreated
     * clinical infection (D) in the next time step.
     */
     if((randNum > phi*fT) && (randNum <= phi)){
       human["state"] = "D";
       human["daysLatent"] = 0;
     }
     /* Asymptomatic infection (E -> A):
      * If the random number is greater than phi, that
      * individual develops an asymptomatic infection (A) in
      * the next time step.
      */
      if(randNum > phi){
        human["state"] = "A";
        human["daysLatent"] = 0;
      }
  }
};

/* T: treated clinical disease */

// [[Rcpp::export]]
void T_compartment(Rcpp::List& human, const Rcpp::NumericVector& theta){

  int dT = Rcpp::as<double>(theta["dT"]);

  double randNum = R::runif(0.0,1.0);
  /* Prophylactic protection (T -> P):
   * If the random number is less than 1/dT, that individual enters the
   * phase of prophylactic protection (P) in the next time step.
  */
  if(randNum <= (1.0/dT)){
    human["state"] = "P";
  }

};

/* D: untreated clinical disease */

// [[Rcpp::export]]
void D_compartment(Rcpp::List& human, const Rcpp::NumericVector& theta){

  int dD = Rcpp::as<double>(theta["dD"]);

  double randNum = R::runif(0.0,1.0);
  /* Progression from diseased to asymptomatic (D -> A):
   * If the random number is less than 1/dD, that individual enters the
   * phase of asymptomatic patent infection (A) in the next time step.
  */
  if(randNum <= (1.0/dD)){
    /* simulation  */
    human["state"] = "A";
  }

};

/* A: asymptomatic patent (detectable by microscopy) infection */

// [[Rcpp::export]]
void A_compartment(Rcpp::List& human, const Rcpp::NumericVector& theta){

  double fT = Rcpp::as<double>(theta["fT"]);
  int dA = Rcpp::as<double>(theta["dA"]);

  double randNum = R::runif(0.0,1.0);

  double phi = Rcpp::as<double>(human["phi"]);
  double lambda = Rcpp::as<double>(human["lambda"]);
  /* Treated clinical infection (A -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
  if(randNum <= phi*fT*lambda){
    human["state"] = "T";
  }
  /* Untreated clinical infection (A -> D):
   * If the random number is greater than phi*fT*lambda and
   * less than phi*lambda, that individual develops an
   * untreated clinical infection (D) in the next time step.
   */
   if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
     human["state"] = "D";
   }
   /* Progression to asymptomatic sub-patent infection (A -> U):
    * If the random number is greater than phi*lambda and less
    * than (phi*lambda + 1/dA), that individual develops sub-patent asymptomatic
    * infection (U) in the next time step.
    */
    if((randNum > phi*lambda) && (randNum <= (phi*lambda + (1.0/dA)))) {
      human["state"] = "U";
    }
};

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */

// [[Rcpp::export]]
void U_compartment(Rcpp::List& human, const Rcpp::NumericVector& theta){

  double fT = Rcpp::as<double>(theta["fT"]);
  int dU = Rcpp::as<double>(theta["dU"]);

  double randNum = R::runif(0.0,1.0);

  double phi = Rcpp::as<double>(human["phi"]);
  double lambda = Rcpp::as<double>(human["lambda"]);

  /* Treated clinical infection (U -> T):
   * If the random number is less than phi*fT*lambda, that
   * individual develops a treated clinical infection (T) in
   * the next time step.
   */
   if(randNum <= phi*fT*lambda){
     human["state"] = "T";
   }
   /* Untreated clinical infection (U -> D):
    * If the random number is greater than phi*fT*lambda and
    * less than phi*lambda, that individual develops an
    * untreated clinical infection (D) in the next time step.
    */
    if((randNum > phi*fT*lambda) && (randNum <= phi*lambda)){
      human["state"] = "D";
    }
    /* Asymptomatic infection (U -> A):
     * If the random number is greater than phi*lambda and
     * less than lambda, that individual develops a patent
     * asymptomatic infection (A) in the next time step.
     */
     if((randNum > phi*lambda) && (randNum <= lambda)){
       human["state"] = "A";
     }
     /* Recovery to susceptible (U -> S):
      * If the random number is greater than lambda and less
      * than (lambda + 1/dU), that individual returns to the susceptible
      * state (S) in the next time step.
      */
      if((randNum > lambda) && (randNum <= (lambda + (1.0/dU)))){
        human["state"] = "S";
      }
};

/* P: protection due to chemoprophylaxis treatment */

// [[Rcpp::export]]
void P_compartment(Rcpp::List& human, const Rcpp::NumericVector& theta){

  int dP = Rcpp::as<double>(theta["dP"]);

  double randNum = R::runif(0.0,1.0);
    /* Prophylactic protection (P -> S):
     * If the random number is less than 1/dP, that individual returns to
     * the susceptible state (S) in the next time step.
     */
    if(randNum <= (1.0/dP)){
      human["state"] = "S";
    }

};


/* ################################################################################
#   Immunity Functions
################################################################################ */


/* immunity */

// [[Rcpp::export]]
void update_immunity(Rcpp::List& human, const Rcpp::NumericVector& theta){

  double uB = Rcpp::as<double>(theta["uB"]);
  double uC = Rcpp::as<double>(theta["uC"]);
  double uD = Rcpp::as<double>(theta["uD"]);
  double dB = Rcpp::as<double>(theta["dB"]);
  double dC = Rcpp::as<double>(theta["dC"]);
  double dID = Rcpp::as<double>(theta["dID"]);
  double dM = Rcpp::as<double>(theta["dM"]);
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
   double epsilon = Rcpp::as<double>(human["epsilon"]);
   double lambda = Rcpp::as<double>(human["lambda"]);

   double IB = Rcpp::as<double>(human["IB"]);
   double ICA = Rcpp::as<double>(human["ICA"]);
   double ICM = Rcpp::as<double>(human["ICM"]);
   double ID = Rcpp::as<double>(human["ID"]);

   human["IB"] = IB + (epsilon/(epsilon*uB + 1.0)) - (IB)*(1.0/dB);
   human["ICA"] = ICA + (lambda/(lambda*uC + 1.0)) - (ICA)*(1.0/dC);
   human["ICM"] = ICM - (ICM)*(1.0/dM);
   human["ID"] = ID + (lambda/(lambda*uD + 1.0)) - (ID)*(1.0/dID);

};

/* lambda */

// [[Rcpp::export]]
void update_lambda(Rcpp::List& human,const Rcpp::NumericVector& theta, const double psi){

  double a0 = Rcpp::as<double>(theta["a0"]);
  double epsilon0 = Rcpp::as<double>(theta["epsilon0"]);
  double b0 = Rcpp::as<double>(theta["b0"]);
  double b1 = Rcpp::as<double>(theta["b1"]);
  double rho = Rcpp::as<double>(theta["rho"]);
  double IB0 = Rcpp::as<double>(theta["IB0"]);
  double kappaB = Rcpp::as<double>(theta["kappaB"]);

  /* Lambda (the force of infection) is calculated for each individual. It
   * varies according to age and biting heterogeneity group:
   */

   double bitingHet = Rcpp::as<double>(human["bitingHet"]);
   double age = Rcpp::as<double>(human["age"]);
   double IB = Rcpp::as<double>(human["IB"]);

   double epsilon = epsilon0 * bitingHet * (1.0 - rho * exp(-age/a0)) * psi;

   human["epsilon"] = epsilon;
   double b = b0*(b1 + ((1.0-b1)/(1.0 + pow((IB/IB0),kappaB))));
   human["lambda"] = epsilon * b;

};

/* phi */

// [[Rcpp::export]]
void update_phi(Rcpp::List& human,const Rcpp::NumericVector& theta){

  double phi0 = Rcpp::as<double>(theta["phi0"]);
  double phi1 = Rcpp::as<double>(theta["phi1"]);
  double IC0 = Rcpp::as<double>(theta["IC0"]);
  double kappaC = Rcpp::as<double>(theta["kappaC"]);

  /* Phi (the probability of acquiring clinical disease upon infection) is
   * also calculated for each individual. It varies according to immune
   * status:
   */

   double ICA = Rcpp::as<double>(human["ICA"]);
   double ICM = Rcpp::as<double>(human["ICM"]);

   human["phi"] = phi0 * (phi1 + ((1.0 - phi1)/(1 + pow(((ICA+ICM)/IC0),kappaC))));

};

/* q (microscopy) */

// [[Rcpp::export]]
void update_q(Rcpp::List& human,const Rcpp::NumericVector& theta){

  double fD0 = Rcpp::as<double>(theta["fD0"]);
  double aD = Rcpp::as<double>(theta["aD"]);
  double gammaD = Rcpp::as<double>(theta["gammaD"]);
  double d1 = Rcpp::as<double>(theta["d1"]);
  double ID0 = Rcpp::as<double>(theta["ID0"]);
  double kappaD = Rcpp::as<double>(theta["kappaD"]);
  double alphaA = Rcpp::as<double>(theta["alphaA"]);
  double alphaU = Rcpp::as<double>(theta["alphaU"]);

 /* q (the probability that an asymptomatic infection is detected by
  * microscopy) is also calculated for each individual, as well as the
  * probability of detection by PCR for asymptomatic infections in states
  * A (patent) and U (subpatent). This also varies according to immune status:
  */
  double ID = Rcpp::as<double>(human["ID"]);
  double age = Rcpp::as<double>(human["age"]);

  double fD = 1.0 - ((1.0 - fD0)/(1 + pow((age/aD),gammaD)));
  double q = d1 + ((1.0 - d1)/(1 + (fD*pow((ID/ID0),kappaD))*fD));

  human["prDetectAMic"] = q;
  human["prDetectAPCR"] = pow(q,alphaA);
  human["prDetectUPCR"] = pow(q,alphaU);

};

// put the functions in a hash table
static const std::unordered_map<std::string, std::function<void(Rcpp::List&,const Rcpp::NumericVector&)> > state_functions  = {
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

// [[Rcpp::export]]
void one_day_human(Rcpp::List& human,const Rcpp::NumericVector& theta, const double psi){

  // mortality
  mortality(human,theta);

  if(Rcpp::as<bool>(human["alive"])){

    // state update
    std::string state = Rcpp::as<std::string>(human["state"]);
    state_functions.at(state)(human,theta);

    // update age
    double age = Rcpp::as<double>(human["age"]);
    human["age"] = age + (1./365.);

    // update immunity
    update_immunity(human,theta);

    // update lambda
    update_lambda(human,theta,psi);

    // update phi
    update_phi(human,theta);

    // update q
    update_q(human,theta);
  }
};
