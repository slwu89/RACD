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
 #  Pre-erythrocytic immunity model
*/

#include <R.h>
#include <math.h>

/*
  Initial parameters

  # age and heterogeneity
  a0: age-dependent biting heterogeneity
  rho age-dependent biting heterogeneity
  zeta: relative biting rate
  psi:

  # immunity reducing probability of infection
  durB: 1/decay rate
  uB: duration in which immunity is not boosted
  b0: probability of infection with no immunity
  b1 maximum relative reduction
  IB0: scale parameter
  kappaB: shape parameter

  # immunity reducing probability of detection
  durD: 1/decay rate
  uD: duration in which immunity is not boosted

  # immunity reducing probabiltiy of clinical disease
  durC: 1/decay rate
  uC: duration in which immunity is not boosted

  # biting
  epsilon0: mean EIR for adults
*/
static double     parms[15];

#define zeta       parms[0]
#define psi        parms[1]
#define a0         parms[2]
#define rho        parms[3]

#define durB       parms[4]
#define uB         parms[5]
#define durD       parms[6]
#define uD         parms[7]
#define durC       parms[8]
#define uC         parms[9]

#define b0         parms[10]
#define b1         parms[11]
#define IB0        parms[12]
#define kappaB     parms[13]
#define epsilon0   parms[14]

/* initializer */
void init_immune(void (* odeparms)(int *, double *))
{
    int N = 15;
    odeparms(&N, parms);
}

/* derivatives */
void derivs_immune(int *neq, double *t, double *y, double *ydot, double *yout, int*ip){

  double time = *t;
  double epsilon = epsilon0*zeta*(1 - rho*exp(-time/a0))*psi; /*EIR at age a*/
  // double b = b0*(b1 + ((1-b1)/(1 + pow((y[0]/IB0),kappaB)))); /*mosquito to human transmission efficiency*/
  double lambda = epsilon*b0*(b1 + ((1.0-b1)/(1.0 + pow((y[0]/IB0),kappaB)))); /*force of infection at age a*/

  ydot[0] = epsilon/(epsilon*uB + 1.0) - y[0]/durB;
  ydot[1] = lambda/(lambda*uD + 1.0) - y[1]/durD;
  ydot[2] = lambda/(lambda*uC + 1.0) - y[2]/durC;

}
