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
void derivs_infection(int *neq, double *t, double *y, double *ydot, double *yout, int*ip){

  double time = *t;
  if(time<0.1){
    printf("epsilon0 %f zeta %f rho %f time %f a0 %f psi %f\n",epsilon0,zeta,rho,time,a0,psi);
  }

  double ICM = initICA20 * exp(-time/dM);
  double epsilon = epsilon0*zeta*(1.0 - rho*exp(-time/a0))*psi; /*EIR at age a*/
  double b = b0*(b1 + ((1-b1)/(1 + pow((y[6]/IB0),kappaB)))); /*mosquito to human transmission efficiency*/
  double lambda = epsilon*b0*(b1 + ((1.0-b1)/(1.0 + pow((y[6]/IB0),kappaB)))); /*force of infection at age a*/
  double phi = phi0*(phi1 + ((1 - phi1)/(1 + pow(((y[7] + ICM)/IC0),kappaC))));

  printf("time: %f ICM: %f epsilon: %f b: %f lambda: %f phi: %f \n",time,ICM,epsilon,b,lambda,phi);

  ydot[0] = -lambda*y[0] + y[5]/dP + y[4]/dU; /* dprS */
  ydot[1] = phi*fT*lambda*(y[0] + y[3] + y[4]) - y[1]/dT; /* dprT */
  ydot[2] = phi*(1 - fT)*lambda*(y[1] + y[3] + y[4]) - y[2]/dD; /* dprD */
  ydot[3] = (1 - phi)*lambda*(y[0] + y[3] + y[4]) + y[2]/dD - lambda*y[3] - y[3]/dA; /* dprA */
  ydot[4] = y[3]/dA - y[4]/dU - lambda*y[4]; /* dprU */
  ydot[5] = y[1]/dT - y[5]/dP; /* dprP */
  ydot[6] = epsilon/(epsilon*uB + 1) - y[6]/dB; /* dIB */
  ydot[7] = lambda/(lambda*uC + 1) - y[7]/dC; /* dICA */

}
