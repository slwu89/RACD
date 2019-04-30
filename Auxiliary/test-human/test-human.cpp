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

#include <math.h>


/* ################################################################################
#   State Transitions
################################################################################ */


/* S: susceptible */

// [[Rcpp::export]]
void S_compartment(Rcpp::List& human){

  double randNum = R::runif(0.0,1.0);
  double lambda = Rcpp::as<double>(human["lambda"]);

  /* Latent infection (S -> E):
   * If the random number is less than lambda, that individual
   * develops a latent infection (E) in the next time step.
   */
  if(randNum <= lambda){
    human["state"] = "E";
    human["daysLatent"] = 1;
  }
};


/* E: latent period */

// [[Rcpp::export]]
void E_compartment(Rcpp::List& human, const int dE, const double fT){
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
void T_compartment(Rcpp::List& human, const int dT){

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
void D_compartment(Rcpp::List& human, const int dD){

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
void A_compartment(Rcpp::List& human, const int dA, const double fT){

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
void U_compartment(Rcpp::List& human, const int dU, const double fT){

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
void P_compartment(Rcpp::List& human, const int dP){

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
void update_immunity(Rcpp::List& human,
  const double uB,
  const double uC,
  const double uD,
  const double dB,
  const double dC,
  const double dID,
  const double dM
){

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
void update_lambda(
  Rcpp::List& human,
  const double a0,
  const double epsilon0,
  const double b0,
  const double b1,
  const double rho,
  const double IB0,
  const double kappaB,
  const double psi
){

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
void update_phi(
  Rcpp::List& human,
  const double phi0,
  const double phi1,
  const double IC0,
  const double kappaC
){

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
void update_q(
  Rcpp::List& human,
  const double fD0,
  const double aD,
  const double gammaD,
  const double d1,
  const double ID0,
  const double kappaD,
  const double alphaA,
  const double alphaU
){

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
