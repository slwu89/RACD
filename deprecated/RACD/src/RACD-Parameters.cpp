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
 #  Parameters Singleton Class Implementation
*/

#include "RACD-Parameters.hpp"

/* constructor & destructor */
parameters::parameters(
  double epsilon0_,
  double fT_,
  int dE_,
  int dT_,
  int dD_,
  int dA_,
  int dU_,
  int dP_,
  double cD_,
  double cT_,
  double cU_,
  double gammaI_,
  double rho_,
  double a0_,
  double sigma2_,
  double d1_,
  double dID_,
  double ID0_,
  double kappaD_,
  double uD_,
  double aD_,
  double fD0_,
  double gammaD_,
  double alphaA_,
  double alphaU_,
  double b0_,
  double b1_,
  double dB_,
  double IB0_,
  double kappaB_,
  double uB_,
  double phi0_,
  double phi1_,
  double dC_,
  double IC0_,
  double kappaC_,
  double uC_,
  double PM_,
  double dM_,
  double rW_,
  double rP_,
  double meanAge_,
  int N_
) : epsilon0(epsilon0_), fT(fT_),
    dE(dE_), dT(dT_), dD(dD_), dA(dA_), dU(dU_), dP(dP_),
    cD(cD_), cT(cT_), cU(cU_), gammaI(gammaI_), rho(rho_),
    a0(a0_), sigma2(sigma2_), d1(d1_), dID(dID_), ID0(ID0_),
    kappaD(kappaD_), uD(uD_), aD(aD_), fD0(fD0_), gammaD(gammaD_),
    alphaA(alphaA_), alphaU(alphaU_), b0(b0_), b1(b1_), dB(dB_),
    IB0(IB0_), kappaB(kappaB_), uB(uB_), phi0(phi0_), phi1(phi1_),
    dC(dC_), IC0(IC0_), kappaC(kappaC_), uC(uC_), PM(PM_),
    dM(dM_), rW(rW_), rP(rP_), meanAge(meanAge_), N(N_), mu(1./(meanAge_*365.))
{
  #ifdef DEBUG_RACD
  std::cout << "parameters being born at " << this << std::endl;
  #endif
};

parameters::~parameters(){
  #ifdef DEBUG_RACD
  std::cout << "parameters being killed at " << this << std::endl;
  #endif
};
