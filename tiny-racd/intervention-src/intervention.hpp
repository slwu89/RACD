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

#ifndef int_mgr_hpp
#define int_mgr_hpp

#include <vector>
#include <memory>

#include <Rcpp.h> // include because we just store a Rcpp::NumericMatrix as the distance matrix btwn houses

// regardless of the type of intervention in play, all the managers need to have the houses
struct house;
using house_ptr = std::unique_ptr<house>;
using house_vector = std::vector<house_ptr>;


/* ################################################################################
#   abstract base intervention manager
################################################################################ */

// the base intervention manager class
class intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_);
  virtual ~intervention_manager() = 0;

  /* factory method (stamp out the type of intervention we want) */
  // NOTE: IF IT DOESNT COMPILE DO THE CONSTRUCTOR INLINE!!!!!!
  static std::unique_ptr<intervention_manager> factory(int type, const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_);

  /* move operators */
  intervention_manager(intervention_manager&&) = default;
  intervention_manager& operator=(intervention_manager&&) = default;

  /* copy operators */
  intervention_manager(intervention_manager&) = delete;
  intervention_manager& operator=(intervention_manager&) = delete;

  /* generic methods */
  virtual void                        one_day_intervention() = 0;

  // after assigning tomorrow's interventions, zero out the data structures tracking them
  void                                zero_house_data();

  // this is for humans to tell the intervention_manager that a clinical case popped up at their house today
  void                                add_cinc(size_t h);

  // history
  Rcpp::IntegerMatrix                 get_int_status_hist(){return int_status_hist;}

protected:

  // track the houses (const pointer to the vec)
  house_vector* const                 houses;
  size_t                              nh;

  // distance matrix between houses
  Rcpp::NumericMatrix                 dmat;

  // radius to search within around each house for interventions
  double                              radius;

  // did each house experience clinical cases today? (ie; does a team need to go out
  // here and "start intervening") on a radius with this house as the center?
  std::vector<bool>                   house_cc;

  // assign intervention (so we don't "intervene" twice; saves computation)
  std::vector<bool>                   house_int;

  // history tracking things
  Rcpp::IntegerMatrix                 int_status_hist;

};


/* ################################################################################
#   RfMDA: reactive focal mass-drug administration
################################################################################ */

class intervention_manager_rfmda : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_rfmda(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_);
  ~intervention_manager_rfmda();

  // implement the RfMDA method
  virtual void one_day_intervention();

private:
};


/* ################################################################################
#   RfVC: reactive focal vector control
################################################################################ */

class intervention_manager_rfvc : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_rfvc(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_);
  ~intervention_manager_rfvc();

  // implement the rfVC method
  virtual void one_day_intervention();

private:
};


/* ################################################################################
#   RACD w/PCR: reactive case detection using PCR
################################################################################ */

class intervention_manager_racd_pcr : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_racd_pcr(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_);
  ~intervention_manager_racd_pcr();

  // implement the RACD w/pcr method
  virtual void one_day_intervention();

private:
};


/* ################################################################################
#   RACD w/Mic: reactive case detection using microscopy
################################################################################ */

class intervention_manager_racd_Mic : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_racd_Mic(const size_t tmax, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_);
  ~intervention_manager_racd_Mic();

  // implement the RACD w/Mic method
  virtual void one_day_intervention();

private:
};


/* ################################################################################
#   RACD w/LAMP: reactive case detection using LAMP
################################################################################ */


#endif