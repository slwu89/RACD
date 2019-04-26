// /*
//  #      ____  ___   __________
//  #     / __ \/   | / ____/ __ \
//  #    / /_/ / /| |/ /   / / / /
//  #   / _, _/ ___ / /___/ /_/ /
//  #  /_/ |_/_/  |_\____/_____/
//  #
//  #  Sean Wu & John M. Marshall
//  #  September 2018
//  #
//  #  Anopheles life cycle model
// */
//
// #include <R.h>
// #include <math.h>
//
//
// /* ################################################################################
//   States:
//   EL: Pre-erythrocytic immunity, reduces the probability of infection
//       following an infectious challenge
//   ID: Detection immunity, blood-stage immunity, reduces the
//       probability of detection and reduces infectiousness to mosquitoes
//   ICA: Acquired clinical immunity, reduces the probability of clinical
//       disease, acquired from previous exposure
//
//   Parameters:
//   a0: age-dependent biting heterogeneity
//   rho age-dependent biting heterogeneity
//   zeta: relative biting rate
//   psi: geographic risk (household level biting heterogeneity)
//
//   # immunity reducing probability of infection
//   durB: 1/decay rate
//   uB: duration in which immunity is not boosted
//   b0: probability of infection with no immunity
//   b1 maximum relative reduction
//   IB0: scale parameter
//   kappaB: shape parameter
//
//   # immunity reducing probability of detection
//   durD: 1/decay rate
//   uD: duration in which immunity is not boosted
//
//   # immunity reducing probabiltiy of clinical disease
//   durC: 1/decay rate
//   uC: duration in which immunity is not boosted
//
//   # biting
//   epsilon0: mean EIR for adults
// ################################################################################ */
//
// static double      parms[15];
//
// #define beta       parms[0]
// #define muEL       parms[1]
// #define muLL       parms[2]
// #define rho        parms[3]
//
// #define durB       parms[4]
// #define uB         parms[5]
// #define durD       parms[6]
// #define uD         parms[7]
// #define durC       parms[8]
// #define uC         parms[9]
//
// #define b0         parms[10]
// #define b1         parms[11]
// #define IB0        parms[12]
// #define kappaB     parms[13]
// #define epsilon0   parms[14]
//
// // beta <- theta[["beta"]] # Eggs laid per day by female mosquito
// // muEL <- theta[["muEL"]] # Early instar stage daily mortality
// // muLL <- theta[["muLL"]] # Late instar stage daily mortality
// // muPL <- theta[["muPL"]] # Pupal stage daily mortality
// // muV <- theta[["muV"]] # Adult mosquito daily mortality
// // Q0 <- theta[["Q0"]] # Human blood index
// // phiB <- theta[["phiB"]] # Proportion of bites on a person while they are in bed
// // phiI <- theta[["phiI"]] # Proportion of bites on a person while they are indoors
// // durEL <- theta[["durEL"]] # Duration of early instar stage (days)
// // durLL <- theta[["durLL"]] # Duration of late instar stage (days)
// // durPL <- theta[["durPL"]] # Duration of pupal stage (days)
// // durEV <- theta[["durEV"]] # Duration of latent period in mosquito (days)
// // gamma <- theta[["gamma"]] # Effect of density-dependence on late instarts relative to early instars
// // tau1 <- theta[["tau1"]] # Time spent foraging for a blood meal (no ITNs) (days)
// // tau2 <- theta[["tau1"]] # Time spent resting and ovipositing (days)
// // NV_eq <- theta[["NV_eq"]] # Number of female mosquitoes at equilibrium
// // lambdaV <- theta[["lambdaV"]] # Force of infection in vectors at equilibrium
// //
// // ## Parameters (interventions):
// // ITNcov <- theta[["ITNcov"]] # ITN coverage
// // IRScov <- theta[["IRScov"]] # IRS coverave
// // time_ITN_on <- theta[["time_ITN_on"]] # When ITNs are applied (days)
// // time_IRS_on <- theta[["time_IRS_on"]] # When IRS is applied (days)
// // rITN <- theta[["rITN"]] # Probability of mosquito repeating a feeding attempt due to IRS
// // sITN <- theta[["sITN"]] # Probability of mosquito feeding and surviving in presence of ITNs
// // rIRS <- theta[["rIRS"]] # Probability of mosquito repeating a feeding attempt due to IRS
// // sIRS <- theta[["sIRS"]] # Probability of mosquito feeding and surviving in presence of IRS
//
//
//
// /* initializer */
// void init_immune(void (* odeparms)(int *, double *))
// {
//     int N = 15;
//     odeparms(&N, parms);
// }
//
// /* derivatives */
// void derivs_immune(int *neq, double *t, double *y, double *ydot, double *yout, int *ip){
//
//   double time = *t;
//
//   /* states */
//   double IB = y[0];
//   double ID = y[1];
//   double ICA = y[2];
//
//   /* parameters */
//   double epsilon = epsilon0*zeta*(1 - rho*exp(-time/a0))*psi; /*EIR at age a*/
//   // double b = b0*(b1 + ((1-b1)/(1 + std::pow((y[0]/IB0),kappaB)))); /*mosquito to human transmission efficiency*/
//   double lambda = epsilon*b0*(b1 + ((1.0-b1)/(1.0 + std::pow((IB/IB0),kappaB)))); /*force of infection at age a*/
//
//   ydot[0] = epsilon/(epsilon*uB + 1.0) - IB/durB;
//   ydot[1] = lambda/(lambda*uD + 1.0) - ID/durD;
//   ydot[2] = lambda/(lambda*uC + 1.0) - ICA/durC;
//
// }
