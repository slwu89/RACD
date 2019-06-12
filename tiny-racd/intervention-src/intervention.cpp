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

// global stuff
#include "globals.hpp"

// other headers we need
#include "house.hpp"
#include "human.hpp"


/* ################################################################################
#   abstract base intervention manager
################################################################################ */

intervention_manager::intervention_manager(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  houses(houses_), nh(nh_), dmat(dmat_), radius(radius_), house_cc(nh), house_int(nh), int_status_hist(tmax,nh) {};

/* define virtual destructor is ok */
intervention_manager::~intervention_manager(){};

// /* move operators */
// intervention_manager::intervention_manager(intervention_manager&&) = default;
// intervention_manager& intervention_manager::operator=(intervention_manager&&) = default;

/* factory method */
std::unique_ptr<intervention_manager> intervention_manager::factory(int type, const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_){

  // 0: RfMDA, 1: RfVC, 2: RACD w/PCR, 3: RACD w/Mic, 4: RACD w/LAMP
  if(type == 0){
    Rcpp::Rcout << "intervention strategy set to: RfMDA\n";
    return std::make_unique<intervention_manager_rfmda>(tmax,houses_,nh_,dmat_,radius_);
  } else if(type == 1){
    Rcpp::Rcout << "intervention strategy set to: RfVC\n";
    return std::make_unique<intervention_manager_rfvc>(tmax,houses_,nh_,dmat_,radius_);
  } else if(type == 2){
    Rcpp::Rcout << "intervention strategy set to: RACD w/PCR\n";
    return std::make_unique<intervention_manager_racd_pcr>(tmax,houses_,nh_,dmat_,radius_);
  } else if(type == 3){
    Rcpp::Rcout << "intervention strategy set to: RACD w/Mic\n";
    return std::make_unique<intervention_manager_racd_Mic>(tmax,houses_,nh_,dmat_,radius_);
  } else if (type == 4){
    Rcpp::stop("RACD w/LAMP not implemented yet!");
  } else{
    Rcpp::stop("intervention type must be an integer in (0,1,2,3,4)");
  }

};

// after assigning tomorrow's interventions, zero out the data structures tracking them
void intervention_manager::zero_house_data(){
  for(size_t h=0; h<nh; h++){
    // house h was a centroid (and by def. was intervened upon)
    if(house_cc[h]){
      int_status_hist.at(tnow,h) = 1;
    }
    // house h was not a centroid but intervened upon because was nearby
    if(!house_cc[h] && house_int[h]){
      int_status_hist.at(tnow,h) = 2;
    }
    // clear out for tomorrow
    house_cc[h] = false;
    house_int[h] = false;
  }
}

// this is for humans to tell the intervention_manager that a clinical case popped up at their house today
void intervention_manager::add_cinc(size_t h){
  house_cc[h] = true;
};


/* ################################################################################
#   RfMDA: reactive focal mass-drug administration
################################################################################ */

/* constructor & destructor */
intervention_manager_rfmda::intervention_manager_rfmda(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(tmax,houses_,nh_,dmat_,radius_) {};

intervention_manager_rfmda::~intervention_manager_rfmda(){};

// implement the RfMDA method
void intervention_manager_rfmda::one_day_intervention(){

  // outermost loop
  for(int h=0; h<nh; h++){

    // if there was no clinical incidence here today, skip it
    if(!house_cc[h]){
      continue;
    }

    // else this house is in the set of centroids for intervention radii
    // first apply the intervention to the people living here if it hasn't been intervened upon already
    if(!house_int[h]){
      apply_MDA(houses->at(h));
      house_int[h] = true;
    }

    // intervene on this house's neighbors
    for(int h_n=0; h_n<nh; h_n++){

      // dont intervene on ourselves or houses too far away or houses which have already been intervened upon
      if((h == h_n) || (dmat.at(h,h_n) > radius) || house_int[h_n]){
        continue;
      }

      // do the intervention
      apply_MDA(houses->at(h_n));
      house_int[h_n] = true;

    }

  }

};


/* ################################################################################
#   RfVC: reactive focal vector control
################################################################################ */

/* constructor & destructor */
intervention_manager_rfvc::intervention_manager_rfvc(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(tmax,houses_,nh_,dmat_,radius_) {};

intervention_manager_rfvc::~intervention_manager_rfvc(){};

// implement the rfVC method
void intervention_manager_rfvc::one_day_intervention(){

  // main loop
  for(int h=0; h<nh; h++){

    // if house h is not in the set of houses where clinical incident cases were
    // picked up today, skip it
    if(!house_cc[h]){
      continue;
    }

    // only intervene on this house if it has not been intervened upon already
    // eg; if it was in the set of neighbors for a previous centroid house
    if(!house_int[h]){
      apply_IRS(houses->at(h));
      house_int[h] = true;
    }

    // household h is the centroid of a circle of radius "radius"; find all other
    // houses in the set of houses in that centroid and intervene upon them
    for(int h_n=0; h_n<nh; h_n++){

      // don't intervene on ourselves, or houses outside the radius, or houses that already got intervention
      if((h == h_n) || (dmat.at(h,h_n) > radius) || house_int[h_n]){
        continue;
      }

      // intervene here
      apply_IRS(houses->at(h_n));
      house_int[h_n] = true;

    }

  }

};


/* ################################################################################
#   RACD w/PCR: reactive case detection using PCR
################################################################################ */

/* constructor & destructor */
intervention_manager_racd_pcr::intervention_manager_racd_pcr(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(tmax,houses_,nh_,dmat_,radius_) {};

intervention_manager_racd_pcr::~intervention_manager_racd_pcr(){};

// implement the RACD w/PCR method
void intervention_manager_racd_pcr::one_day_intervention(){

  // main loop
  for(int h=0; h<nh; h++){

    // if house h is not in the set of houses where clinical incident cases were
    // picked up today, skip it
    if(!house_cc[h]){
      continue;
    }

    // only intervene on this house if it has not been intervened upon already
    // eg; if it was in the set of neighbors for a previous centroid house
    if(!house_int[h]){
      apply_RACD_PCR(houses->at(h));
      house_int[h] = true;
    }

    // household h is the centroid of a circle of radius "radius"; find all other
    // houses in the set of houses in that centroid and intervene upon them
    for(int h_n=0; h_n<nh; h_n++){

      // don't intervene on ourselves, or houses outside the radius, or houses that already got intervention
      if((h == h_n) || (dmat.at(h,h_n) > radius) || house_int[h_n]){
        continue;
      }

      // intervene here
      apply_RACD_PCR(houses->at(h_n));
      house_int[h_n] = true;

    }

  }

};


/* ################################################################################
#   RACD w/Mic: reactive case detection using microscopy
################################################################################ */

/* constructor & destructor */
intervention_manager_racd_Mic::intervention_manager_racd_Mic(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_) :
  intervention_manager(tmax,houses_,nh_,dmat_,radius_) {};

intervention_manager_racd_Mic::~intervention_manager_racd_Mic(){};

// implement the RACD w/Mic method
void intervention_manager_racd_Mic::one_day_intervention(){

  // main loop
  for(int h=0; h<nh; h++){

    // if house h is not in the set of houses where clinical incident cases were
    // picked up today, skip it
    if(!house_cc[h]){
      continue;
    }

    // only intervene on this house if it has not been intervened upon already
    // eg; if it was in the set of neighbors for a previous centroid house
    if(!house_int[h]){
      apply_RACD_Mic(houses->at(h));
      house_int[h] = true;
    }

    // household h is the centroid of a circle of radius "radius"; find all other
    // houses in the set of houses in that centroid and intervene upon them
    for(int h_n=0; h_n<nh; h_n++){

      // don't intervene on ourselves, or houses outside the radius, or houses that already got intervention
      if((h == h_n) || (dmat.at(h,h_n) > radius) || house_int[h_n]){
        continue;
      }

      // intervene here
      apply_RACD_Mic(houses->at(h_n));
      house_int[h_n] = true;

    }

  }

};


/* ################################################################################
#   RACD w/LAMP: reactive case detection using LAMP
################################################################################ */
