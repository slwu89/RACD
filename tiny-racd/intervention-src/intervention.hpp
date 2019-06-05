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
  house_vector* const       houses;

};


/* ################################################################################
#   RfMDA: reactive focal mass-drug administration
################################################################################ */


/* ################################################################################
#   RfVC: reactive focal vector control
################################################################################ */


/* ################################################################################
#   RACD w/PCR: reactive case detection using PCR
################################################################################ */


/* ################################################################################
#   RACD w/Mic: reactive case detection using microscopy
################################################################################ */


/* ################################################################################
#   RACD w/LAMP: reactive case detection using LAMP
################################################################################ */


#endif
