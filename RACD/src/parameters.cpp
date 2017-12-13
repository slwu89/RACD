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

#include "parameters.hpp"
#include "prng.hpp"

/* declare pointer to instance */
RACD_Parameters* RACD_Parameters::RACD_Parameters_instance = nullptr;

/* constructor */
RACD_Parameters::RACD_Parameters(){
  #ifdef DEBUG_HPP
  std::cout << "RACD_Parameters being born at " << this << std::endl;
  #endif
};

/* destructor */
RACD_Parameters::~RACD_Parameters(){
  #ifdef DEBUG_HPP
  std::cout << "RACD_Parameters being killed at " << this << std::endl;
  #endif
};

/* set parameter values */
void RACD_Parameters::set_values(
  double _epsilon0,
  double _fT,
  int _dE,
  int _dT,
  int _dD,
  int _dA,
  int _dU,
  int _dP,
  double _cD,
  double _cT,
  double _cU,
  double _gammaI,
  double _rho,
  double _a0,
  double _sigma2,
  double _d1,
  double _dID,
  double _ID0,
  double _kappaD,
  double _uD,
  double _aD,
  double _fD0,
  double _gammaD,
  double _alphaA,
  double _alphaU,
  double _b0,
  double _b1,
  double _dB,
  double _IB0,
  double _kappaB,
  double _uB,
  double _phi0,
  double _phi1,
  double _dC,
  double _IC0,
  double _kappaC,
  double _uC,
  double _PM,
  double _dM,
  double _rW,
  double _rP,
  double _meanAge,
  int _N,
  double _meanNumPeoplePerHouse,
  int _numHousesPerBreedingSite
){
  epsilon0 = _epsilon0;
  fT = _fT;
  dE = _dE;
  dT = _dT;
  dD = _dD;
  dA = _dA;
  dU = _dU;
  dP = _dP;
  cD = _cD;
  cT = _cT;
  cU = _cU;
  gammaI = _gammaI;
  rho = _rho;
  a0 = _a0;
  sigma2 = _sigma2;
  d1 = _d1;
  dID = _dID;
  ID0 = _ID0;
  kappaD = _kappaD;
  uD = _uD;
  aD = _aD;
  fD0 = _fD0;
  gammaD = _gammaD;
  alphaA = _alphaA;
  alphaU = _alphaU;
  b0 = _b0;
  b1 = _b1;
  dB = _dB;
  IB0 = _IB0;
  kappaB = _kappaB;
  uB = _uB;
  phi0 = _phi0;
  phi1 = _phi1;
  dC = _dC;
  IC0 = _IC0;
  kappaC = _kappaC;
  uC = _uC;
  PM = _PM;
  dM = _dM;
  rW = _rW;
  rP = _rP;
  meanAge = _meanAge;
  N = _N;
  meanNumPeoplePerHouse = _meanNumPeoplePerHouse;
  numHousesPerBreedingSite = _numHousesPerBreedingSite;
};

/* pseudo-random number generation */
void RACD_Parameters::set_prng(const uint_least32_t &seed) {
  prng_member = std::make_unique<prng>(seed);
};

prng* RACD_Parameters::get_prng(){
  return prng_member.get();
};

/* suicide at end of program */
void RACD_Parameters::suicide(){
  #ifdef DEBUG_HPP
  std::cout << "RACD_Parameters being cleared at " << this << std::endl;
  #endif
  RACD_Parameters::~RACD_Parameters();
};

/* return instance */
RACD_Parameters* RACD_Parameters::instance(){
  if (!RACD_Parameters_instance)
    RACD_Parameters_instance = new RACD_Parameters;
  return RACD_Parameters_instance;
};
