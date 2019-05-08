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
 #  the human
*/

#ifndef human_hpp
#define human_hpp

#include <iostream>

// we're all special snowflakes here
static size_t global_hid = 0;

// a person
typedef struct human {

  // known on initialization
  size_t        id;
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
  double        c;
  std::string   state;

  // other dynamic members
  size_t         days_latent;
  bool           ITN;
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

// the smart pointer for a human
using human_ptr = std::unique_ptr<human>;










#endif
