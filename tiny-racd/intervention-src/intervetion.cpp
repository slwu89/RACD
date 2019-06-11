/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Marshall Lab (https://www.marshalllab.com)
 #  Sean Wu (slwu89@berkeley.edu)
 #  May 2019
 #
 #  intervention manager class
*/

#include "intervention.hpp"

// other headers we need
#include "house.hpp"
#include "human.hpp"


/* ################################################################################
#   abstract base intervention manager
################################################################################ */

intervention_manager::intervention_manager(house_vector* houses_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  houses(houses_), dmat(dmat_), radius(radius_), house_cc(houses_->size()), house_int(houses_->size()) {};

/* define virtual destructor is ok */
intervention_manager::~intervention_manager(){};

// /* move operators */
// intervention_manager::intervention_manager(intervention_manager&&) = default;
// intervention_manager& intervention_manager::operator=(intervention_manager&&) = default;

/* factory method */
std::unique_ptr<intervention_manager> factory(int type){

  // /* check what model we need to make */
  // std::string model(Rcpp::as<std::string>(human_pars["model"]));
  //
  // /* make a derived class */
  // if(model.compare("PfSI") == 0){
  //   return std::make_unique<human_pfsi>(human_pars,tileP_);
  // } else {
  //   Rcpp::stop("invalid 'model' field in 'human_pars'\n");
  // }
};



/* ################################################################################
#   RfMDA: reactive focal mass-drug administration
################################################################################ */

/* constructor & destructor */
intervention_manager_rfmda::intervention_manager_rfmda(house_vector* houses_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(houses_,dmat_,radius_) {};

intervention_manager_rfmda::~intervention_manager_rfmda(){};

  // implement the rfVC method
void intervention_manager_rfmda::one_day_intervention(){

};


/* ################################################################################
#   RfVC: reactive focal vector control
################################################################################ */

/* constructor & destructor */
intervention_manager_rfvc::intervention_manager_rfvc(house_vector* houses_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(houses_,dmat_,radius_) {};

intervention_manager_rfvc::~intervention_manager_rfvc(){};


/* ################################################################################
#   RACD w/PCR: reactive case detection using PCR
################################################################################ */

/* constructor & destructor */
intervention_manager_racd_pcr::intervention_manager_racd_pcr(house_vector* houses_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(houses_,dmat_,radius_) {};

intervention_manager_racd_pcr::~intervention_manager_racd_pcr(){};


/* ################################################################################
#   RACD w/Mic: reactive case detection using microscopy
################################################################################ */

/* constructor & destructor */
intervention_manager_racd_Mic::intervention_manager_racd_Mic(house_vector* houses_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(houses_,dmat_,radius_) {};

intervention_manager_racd_Mic::~intervention_manager_racd_Mic(){};


/* ################################################################################
#   RACD w/LAMP: reactive case detection using LAMP
################################################################################ */
