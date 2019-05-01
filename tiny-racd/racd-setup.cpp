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
 #  setup RACD
*/

#include <Rcpp.h>
#include <math.h>

// EIR_h: my personal EIR (EIR_tot * pi)

// [[Rcpp::export]]
Rcpp::List immmune_ode(
  const double time,
  const Rcpp::NumericVector& state,
  const Rcpp::NumericVector& theta,
  const double EIR_h
){
  // states
  double IB = state[0];
  double ID = state[1];
  double ICA = state[2];

  // dx/dt
  Rcpp::NumericVector dx(3);

  // params
  double durB = Rcpp::as<double>(theta["dB"]) / 365.; // Inverse decay rate (years)
  double uB = Rcpp::as<double>(theta["uB"]) / 365.; // Duration in which immunity is not boosted (years)
  double durD = Rcpp::as<double>(theta["dID"]) / 365.; // Inverse decay rate (years)
  double uD = Rcpp::as<double>(theta["uD"]) / 365.; // Duration in which immunity is not boosted (years)
  double durC = Rcpp::as<double>(theta["dC"]) / 365.; // Inverse decay rate (years)
  double uC = Rcpp::as<double>(theta["uC"]) / 365.; // Duration in which immunity is not boosted (years)
  double b0 = Rcpp::as<double>(theta["b0"]);
  double b1 = Rcpp::as<double>(theta["b1"]);
  double IB0 = Rcpp::as<double>(theta["IB0"]);
  double kappaB = Rcpp::as<double>(theta["kappaB"]);

  double EIR = EIR_h * 365.; // my EIR (years)
  double b = b0*(b1 + ((1.-b1)/(1. + std::pow((IB/IB0),kappaB))));
  double lambda = EIR * b;

  // ODEs
  dx[0] = EIR/(EIR*uB + 1.) - IB/durB; // IB
  dx[1] = lambda/(lambda*uD + 1.) - ID/durD;   // ID
  dx[2] = lambda/(lambda*uC + 1.) - ICA/durC;  // ICA

  return Rcpp::List::create(dx);
};


// EIR_h: my personal EIR (EIR_tot * pi)

// [[Rcpp::export]]
Rcpp::List state_ode(
  const double time,
  const Rcpp::NumericVector& state,
  const Rcpp::NumericVector& theta,
  const double EIR_h
){

  // states
  double prS = state[0];
  double prT = state[1];
  double prD = state[2];
  double prA = state[3];
  double prU = state[4];
  double prP = state[5];
  double IB = state[6];
  double ICA = state[7];

  // dx/dt
  Rcpp::NumericVector dx(8);

  // Parameters:
  double fT = Rcpp::as<double>(theta["fT"]); // Proportion of clinical disease cases successfully treated

  // Model parameters taken from Griffin et al. (2014):
  // Human infection durations:
  double dT = Rcpp::as<double>(theta["dT"]) / 365.; // Duration of treated clinical disease (years)
  double dD = Rcpp::as<double>(theta["dD"]) / 365.; // Duration of untreated clinical disease (years)
  double dA = Rcpp::as<double>(theta["dA"]) / 365.; // Duration of patent infection (years)
  double dU = Rcpp::as<double>(theta["dU"]) / 365.; // Duration of sub-patent infection (years) (fitted)
  double dP = Rcpp::as<double>(theta["dP"]) / 365.; // Duration of prophylactic protection following treatment (years)

  // Immunity reducing probability of infection:
  double b0 = Rcpp::as<double>(theta["b0"]); // Probabiliy with no immunity (fitted)
  double b1 = Rcpp::as<double>(theta["b1"]); // Maximum relative reduction
  double dB = Rcpp::as<double>(theta["dB"]) / 365.; // Inverse of decay rate (years)
  double IB0 = Rcpp::as<double>(theta["IB0"]); // Scale parameter (fitted)
  double kappaB = Rcpp::as<double>(theta["kappaB"]); // Shape parameter (fitted)
  double uB = Rcpp::as<double>(theta["uB"]) / 365.; // Duration in which immunity is not boosted (years) (fitted)

  // Immunity reducing probability of clinical disease:
  double phi0 = Rcpp::as<double>(theta["phi0"]); // Probability with no immunity
  double phi1 = Rcpp::as<double>(theta["phi1"]); // Maximum relative reduction
  double dC = Rcpp::as<double>(theta["dC"]) / 365.; // Inverse decay rate (years)
  double IC0 = Rcpp::as<double>(theta["IC0"]); // Scale parameter
  double kappaC = Rcpp::as<double>(theta["kappaC"]); // Shape parameter
  double uC = Rcpp::as<double>(theta["uC"]) / 365.; // Duration in which immunity is not boosted (years)
  double dM = Rcpp::as<double>(theta["dM"]) / 365.; // Inverse decay rate of maternal immunity (years)
  double initICA20 = Rcpp::as<double>(theta["initICA20"]);

  double ICM = initICA20 * std::exp(-time/dM);

  double EIR = EIR_h * 365.; // my EIR (years)
  double b = b0*(b1 + ((1.-b1)/(1. + std::pow((IB/IB0),kappaB))));
  double lambda = EIR * b;
  double phi = phi0 * (phi1 + ((1. - phi1)/(1. + std::pow(((ICA+ICM)/IC0),kappaC))));

  // ODEs
	dx[0] = - lambda*prS + prP/dP + prU/dU;
	dx[1] = phi*fT*lambda*(prS + prA + prU) - prT/dT;
	dx[2] = phi*(1. - fT)*lambda*(prS + prA + prU) - prD/dD;
	dx[3] = (1. - phi)*lambda*(prS + prA + prU) + prD/dD - lambda*prA - prA/dA;
	dx[4] = prA/dA - prU/dU - lambda*prU;
	dx[5] = prT/dT - prP/dP;
	dx[6] = EIR/(EIR*uB + 1.) - IB/dB;
	dx[7] = lambda/(lambda*uC + 1.) - ICA/dC;

  return Rcpp::List::create(dx);
};
