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

#ifndef HUMAN_HPP
#define HUMAN_HPP

#include <iostream>
#include <memory>
#include <unordered_map>

#include <Rcpp.h>

// forward declare house
struct house;


/* ################################################################################
#   the human
################################################################################ */

// a person
typedef struct human {

  static int    global_id;

  // known on initialization
  int           id;
  double        age;
  bool          alive;
  house*        house_ptr;
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
  int           days_latent;

  // intervention
  bool           ITN;
  int            ITN_time_off;

  // constructor/destructor
  human(const double age_,
        house* house_ptr_,
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
        const double c_,
        const std::string state_);
  ~human();
} human;

int human::global_id = 0;

// the smart pointer for a human
using human_ptr = std::unique_ptr<human>;




/* ################################################################################
#   track output
################################################################################ */

Rcpp::List human_2list(human const* const human_ptr);

#endif
