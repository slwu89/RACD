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

#include <Rcpp.h>

#include <iostream>
#include <memory>
#include <unordered_map>


// forward declaration
struct house;

// a person
typedef struct human {

  // known on initialization
  size_t        id;
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
  size_t         days_latent;

  // intervention
  bool           ITN;
  size_t         ITN_time_off;

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

// the smart pointer for a human
using human_ptr = std::unique_ptr<human>;


/* ################################################################################
#   individual level tracking
################################################################################ */

// track clinical incidence
void track_cinc(const human_ptr& h);


/* ################################################################################
#   State transitions for our little Markov humans
################################################################################ */

void mortality(human_ptr& human);

/* S: susceptible */
void S_compartment(human_ptr& human);

/* E: latent period */
void E_compartment(human_ptr& human);

/* T: treated clinical disease */
void T_compartment(human_ptr& human);

/* D: untreated clinical disease */
void D_compartment(human_ptr& human);

/* A: asymptomatic patent (detectable by microscopy) infection */
void A_compartment(human_ptr& human);

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void U_compartment(human_ptr& human);

/* P: protection due to chemoprophylaxis treatment */
void P_compartment(human_ptr& human);


/* ################################################################################
#   Immunity Functions
################################################################################ */

/* immunity */
void update_immunity(human_ptr& human);

/* lambda */
void update_lambda(human_ptr& human);

/* phi */
void update_phi(human_ptr& human);

/* q (microscopy) */
void update_q(human_ptr& human);

// put the functions in a hash table
static std::unordered_map<std::string, std::function<void(human_ptr& human)> > state_functions  = {
  {"S",S_compartment},
  {"E",E_compartment},
  {"T",T_compartment},
  {"D",D_compartment},
  {"A",A_compartment},
  {"U",U_compartment},
  {"P",P_compartment}
};


/* ################################################################################
#   Infectiousness to Mosquitoes
################################################################################ */

void infectiousness_S(human_ptr& human);

void infectiousness_E(human_ptr& human);

void infectiousness_T(human_ptr& human);

void infectiousness_D(human_ptr& human);

void infectiousness_A(human_ptr& human);

void infectiousness_U(human_ptr& human);

void infectiousness_P(human_ptr& human);

// put the functions in a hash table
static std::unordered_map<std::string, std::function<void(human_ptr& human)> > infectivity_functions  = {
  {"S",infectiousness_S},
  {"E",infectiousness_E},
  {"T",infectiousness_T},
  {"D",infectiousness_D},
  {"A",infectiousness_A},
  {"U",infectiousness_U},
  {"P",infectiousness_P}
};


/* ################################################################################
#   Mosquito Approch probabilities (what happens when a bloodsucker tries to bite me)
################################################################################ */

/* probability of successful biting */
double get_w(human_ptr& human);

/* probability of biting */
double get_y(human_ptr& human);

/* probability of repellency*/
double get_z(human_ptr& human);


/* ################################################################################
#   bookkeeping (bites and biting weights)
#   these work on humans and houses
################################################################################ */

// add my biting weight to the house
void add_pi(human_ptr& human);

// take out my biting weight from the house
void remove_pi(human_ptr& human);

// update pi (do this daily because I age)
void update_pi(human_ptr& human);


/* ################################################################################
#   Humans: daily update
################################################################################ */

void one_day_update_human(human_ptr& human);


/* ################################################################################
#   interventions
################################################################################ */

// called before exiting daily update; check if interventions expire
void update_interventions_human(human_ptr& human);

void give_ITN(human_ptr& human);


#endif
