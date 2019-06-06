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
  intervention_manager(house_vector* houses_);
  virtual ~intervention_manager() = 0;

  /* factory method (stamp out the type of intervention we want) */
  // NOTE: IF IT DOESNT COMPILE DO THE CONSTRUCTOR INLINE!!!!!!
  static std::unique_ptr<intervention_manager> factory(int type);

  /* move operators */
  intervention_manager(intervention_manager&&) = default;
  intervention_manager& operator=(intervention_manager&&) = default;

  /* copy operators */
  intervention_manager(intervention_manager&) = delete;
  intervention_manager& operator=(intervention_manager&) = delete;

protected:

  // track the houses (const pointer to the vec)
  house_vector* const                 houses;

  // distance matrix between houses
  std::vector<std::vector<double> >   dmat;

};


/* ################################################################################
#   RfMDA: reactive focal mass-drug administration
################################################################################ */

class intervention_manager_rfmda : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_rfmda(house_vector* houses_);
  ~intervention_manager_rfmda();

private:
};


/* ################################################################################
#   RfVC: reactive focal vector control
################################################################################ */

class intervention_manager_rfvc : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_rfvc(house_vector* houses_);
  ~intervention_manager_rfvc();

private:
};


/* ################################################################################
#   RACD w/PCR: reactive case detection using PCR
################################################################################ */

class intervention_manager_racd_pcr : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_racd_pcr(house_vector* houses_);
  ~intervention_manager_racd_pcr();

private:
};


/* ################################################################################
#   RACD w/Mic: reactive case detection using microscopy
################################################################################ */

class intervention_manager_racd_Mic : public intervention_manager {
public:

  /* constructor & destructor */
  intervention_manager_racd_Mic(house_vector* houses_);
  ~intervention_manager_racd_Mic();

private:
};


/* ################################################################################
#   RACD w/LAMP: reactive case detection using LAMP
################################################################################ */


#endif
