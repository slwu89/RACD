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
 #  Infection status model
*/

#include <R.h>
#include <math.h>

/* ################################################################################

  Differential equations for calculating the probability that an individual
  is in each state given their age and EIR heterogeneity attributes:

  Parameters:
  Variable model parameters:
  epsilon0: Mean EIR for adults (per year)
  fT: proportion of clinical disease cases successfully treated

  Model parameters taken from Griffin et al. (2014):
  Human infection durations:
  dE: Duration of latent period (years)
  dT: Duration of treated clinical disease (years)
  dD: Duration of untreated clinical disease (years)
  dA: Duration of patent infection (years)
  dU: Duration of sub-patent infection (years) (fitted)
  dP: Duration of prophylactic protection following treatment (years)

  Age parameters:
  rho: Age-dependent biting parameter
  a0: Age-dependent biting parameter (years)

  Immunity reducing probability of infection:
  b0: Probabiliy with no immunity (fitted)
  b1: Maximum relative reduction
  dB: Inverse of decay rate (years)
  IB0: Scale parameter (fitted)
  kappaB: Shape parameter (fitted)
  uB: Duration in which immunity is not boosted (years) (fitted)

  Immunity reducing probability of clinical disease:
  phi0: Probability with no immunity
  phi1: Maximum relative reduction
  dC: Inverse decay rate (years)
  IC0: Scale parameter
  kappaC: Shape parameter
  uC: Duration in which immunity is not boosted (years)
  PM: New-born immunity relative to mother's immunity
  dM: Inverse decay rate of maternal immunity (years)

  initICA20: maternally-derived antibodies

  Individual heterogeneity parameters:
  zeta: individual level biting heterogeneity
  psi: geographic risk (household level biting heterogeneity)

################################################################################ */

static double         parms[27];

/* model parameters */
#define epsilon0      parms[0]
#define fT            parms[1]

/* human infection durations */
#define dE            parms[2]
#define dT            parms[3]
#define dD            parms[4]
#define dA            parms[5]
#define dU            parms[6]
#define dP            parms[7]

/* biting heterogeneity */
#define rho           parms[8]
#define a0            parms[9]

/* immunity reducing probability of infection */
#define b0            parms[10]
#define b1            parms[11]
#define dB            parms[12]
#define IB0           parms[13]
#define kappaB        parms[14]
#define uB            parms[15]

/* immunity reducing probability of clinical disease */
#define phi0          parms[16]
#define phi1          parms[17]
#define dC            parms[18]
#define IC0           parms[19]
#define kappaC        parms[20]
#define uC            parms[21]
#define PM            parms[22]
#define dM            parms[23]
#define initICA20     parms[24]

/* individual heterogeneity parameters */
#define zeta          parms[25]
#define psi           parms[26]

/* initializer */
void init_infection(void (* odeparms)(int *, double *)){
    int N = 27;
    odeparms(&N, parms);
}

/* derivatives */
void derivs_infection(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){

  double time = *t;

  /* states */
  double prS = y[0];
  double prT = y[1];
  double prD = y[2];
  double prA = y[3];
  double prU = y[4];
  double prP = y[5];
  double IB = y[6];
  double ICA = y[7];

  /* parameters */
  double ICM = initICA20 * exp(-time/dM); /*maternally-inherited antibodies*/
  double epsilon = epsilon0*zeta*(1 - rho*exp(-time/a0))*psi; /*EIR at age a*/
  // double b = b0*(b1 + ((1-b1)/(1 + std::pow((y[6]/IB0),kappaB)))); /*mosquito to human transmission efficiency*/
  double lambda = epsilon*b0*(b1 + ((1-b1)/(1 + pow((IB/IB0),kappaB)))); /*force of infection at age a*/
  double phi = phi0*(phi1 + ((1 - phi1)/(1 + pow(((ICA + ICM)/IC0),kappaC)))); /* probability to acquire clinical disease if infected */

  ydot[0] = -lambda*prS + prP/dP + prU/dU; /* dprS */
  ydot[1] = phi*fT*lambda*(prS + prA + prU) - prT/dT; /* dprT */
  ydot[2] = phi*(1 - fT)*lambda*(prS + prA + prU) - prD/dD; /* dprD */
  ydot[3] = (1 - phi)*lambda*(prS + prA + prU) + prD/dD - lambda*prA - prA/dA; /* dprA */
  ydot[4] = prA/dA - prU/dU - lambda*prU; /* dprU */
  ydot[5] = prT/dT - prP/dP; /* dprP */
  ydot[6] = epsilon/(epsilon*uB + 1) - IB/dB; /* dIB */
  ydot[7] = lambda/(lambda*uC + 1) - ICA/dC; /* dICA */

}
