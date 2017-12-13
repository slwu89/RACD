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
 #  Human Class Implementation
*/

#include "RACD-Human.hpp"
#include "RACD-House.hpp"
#include "RACD-Parameters.hpp"
#include "RACD-PRNG.hpp"

#include <vector>

/* constructor */
human::human(const int& _humanID,
      const double& _age,
      const bool& _alive,
      house* _house_ptr,
      const double& _bitingHet,
      const double& _IB,
      const double& _ID,
      const double& _ICA,
      const double& _ICM,
      const double& _epsilon,
      const double& _lambda,
      const double& _phi,
      const double& _prDetectAMic,
      const double& _prDetectAPCR,
      const double& _prDetectUPCR,
      const std::string _state,
      const int _daysLatent
    ) : humanID(_humanID), age(_age), alive(_alive), house_ptr(_house_ptr),
        bitingHet(_bitingHet), IB(_IB), ID(_ID), ICA(_ICA), ICM(_ICM),
        epsilon(_epsilon), lambda(_lambda), phi(_phi),
        prDetectAMic(_prDetectAMic), prDetectAPCR(_prDetectAPCR), prDetectUPCR(_prDetectUPCR),
        state(_state), daysLatent(_daysLatent)
{
  #ifdef DEBUG_HPP
  std::cout << "human " << humanID << " being born at " << this << std::endl;
  #endif
};

/* destructor */
human::~human(){
  #ifdef DEBUG_HPP
  std::cout << "human " << humanID << " being killed at " << this << std::endl;
  #endif
};

/* object suicide */
void human::suicide(){
  #ifdef DEBUG_HPP
  std::cout << "human " << humanID << " suiciding at " << this << std::endl;
  #endif
  for(auto it = house_ptr->get_humans().begin(); it != house_ptr->get_humans().end(); it++){
    if(it->get() == this){
      house_ptr->get_humans().erase(it);
    }
  }
}


/* Simulation Methods */

/* mortality */
void human::mortality(){

  /* Determine if I die during the current time step */
  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();
  if(randNum < RACD_Parameters::instance()->get_mu()){
    // REPLACE WITH REAL OUTPUT
    std::cout << "human " << humanID << " died due to natural mortality" << std::endl;
    // REPLACE WITH REAL OUTPUT
    alive = false;
    suicide();
  }

};

/* S: susceptible */
void human::S_compartment(){

  /* Draw a random number between 1 and 0 for each susceptible individual */
  double randNum = RACD_Parameters::instance()->get_prng()->get_runif();

  /* Latent infection (S -> E):
   * If the random number is less than lambda, that individual
   * develops a latent infection (E) in the next time step.
   */ 
  if(randNum <= lambda){
    state = "E";
    daysLatent = 1;
  }
};


/* E: latent period */
void                            E_compartment();

/* T: treated clinical disease */
void                            T_compartment();

/* D: untreated clinical disease */
void                            D_compartment();

/* A: asymptomatic patent (detectable by microscopy) infection */
void                            A_compartment();

/* U: asymptomatic sub-patent (not detectable by microscopy) infection */
void                            U_compartment();

/* P: protection due to chemoprophylaxis treatment */
void                            P_compartment();
