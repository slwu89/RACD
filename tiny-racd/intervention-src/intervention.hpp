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
  intervention_manager(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int tdelay_);
  virtual ~intervention_manager() = 0;

  /* factory method (stamp out the type of intervention we want) */
  // NOTE: IF IT DOESNT COMPILE DO THE CONSTRUCTOR INLINE!!!!!!
  static std::unique_ptr<intervention_manager> factory(int type, const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int tdelay_);

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

  // time to start intervention
  int                                 tstart;
  int                                 tend;

  // delay in consecutive sweeps
  std::vector<int>                    tlast;
  int                                 tdelay;

  // did each house experience clinical cases today? (ie; does a team need to go out
  // here and "start intervening") on a radius with this house as the center?
  std::vector<bool>                   house_cc;

  // assign intervention (so we don't "intervene" twice; saves computation)
  std::vector<bool>                   house_int;

  // history tracking things
  Rcpp::IntegerMatrix                 int_status_hist;

};


/* ################################################################################
#   NULL: no intervention
################################################################################ */

extern Rcpp::NumericMatrix null_dmat;

class intervention_manager_null : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_null(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int tdelay_);
  ~intervention_manager_null();

  // implement the null method
  virtual void one_day_intervention();

private:
};


/* ################################################################################
#   RfMDA: reactive focal mass-drug administration
################################################################################ */

class intervention_manager_rfmda : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_rfmda(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int tdelay_);
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
  intervention_manager_rfvc(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int max_house_, const int tdelay_);
  ~intervention_manager_rfvc();

  // implement the rfVC method
  virtual void one_day_intervention();

private:
  // the RfVC protocol is to only spray 7 houses maximum in the 500m, so need to stop afterwards
  int   max_house; // max houses to spray
  int   house_count; // current houses sprayed
};


/* ################################################################################
#   RACD w/PCR: reactive case detection using PCR
################################################################################ */

class intervention_manager_racd_pcr : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_racd_pcr(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int tdelay_);
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
  intervention_manager_racd_Mic(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int tdelay_);
  ~intervention_manager_racd_Mic();

  // implement the RACD w/Mic method
  virtual void one_day_intervention();

private:
};


/* ################################################################################
#   RACD w/LAMP: reactive case detection using LAMP
################################################################################ */


/* ################################################################################
#   RACD w/Mic + RfVC:
#   reactive case detection using LAMP + reactive focal vector control
################################################################################ */

class intervention_manager_racdMic_rfvc : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_racdMic_rfvc(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int max_house_, const int tdelay_);
  ~intervention_manager_racdMic_rfvc();

  // implement the rfVC method
  virtual void one_day_intervention();

private:
  // the RfVC protocol is to only spray 7 houses maximum in the 500m, so need to stop afterwards
  int   max_house; // max houses to spray
  int   house_count; // current houses sprayed
};


/* ################################################################################
#   RfMDA + RfVC:
#   reactive focal mass-drug administration + reactive focal vector control
################################################################################ */

class intervention_manager_rfmda_rfvc : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_rfmda_rfvc(const size_t tmax_, const int tstart_, const int tend_, house_vector* houses_, const size_t nh_, const Rcpp::NumericMatrix& dmat_, const double radius_, const int max_house_, const int tdelay_);
  ~intervention_manager_rfmda_rfvc();

  // implement the rfVC method
  virtual void one_day_intervention();

private:
  // the RfVC protocol is to only spray 7 houses maximum in the 500m, so need to stop afterwards
  int   max_house; // max houses to spray
  int   house_count; // current houses sprayed
};



#endif
