/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Marshall Lab
 #  May 2018
 #
 #  Mosquito Habitat Class Definition
*/

#ifndef RACD_HABITAT_
#define RACD_HABITAT_

#include <iostream>
#include <math.h>
#include <memory>

/* alias and forward declarations */
class village;

/* ################################################################################
 * abstract base class
 *  basically, the mosquito habitats need to simulate on a one day time step
 *  and send out bites to houses. They have a pointer to the containing village
 *  that allows them to do this.
 *
################################################################################ */

class mosquito_habitat_base {

public:
  /* define the interface */

  /* constructor */
  mosquito_habitat_base(const int habitatID_, village* village_ptr_);

  /* pure virtual destructor */
  virtual ~mosquito_habitat_base() = 0;

  /* required implementations */
  virtual void one_day(const int tNow) = 0;

  /* accessors */
  int get_habitatID(){ return habitatID; };
  village* get_village_ptr(){ return village_ptr; };

protected:
  /* default data members */
  int                         habitatID;
  village*                    village_ptr; /* raw pointer ok because house lifespan > human lifespan in memory */
};


/* ################################################################################
 * default mosquito model
 *  seasonally forced EIR
################################################################################ */

class mosquito_habitat_eir : public mosquito_habitat_base {

public:

  /* constructor */
  mosquito_habitat_eir(const double EIR_mean_, const double offset_, const int habitatID_, village* village_ptr_);

  /* destructor */
  ~mosquito_habitat_eir();

  /* daily simulation */
  void one_day(const int tNow);
private:
  double                      EIR_mean;
  double                      offset; /**/
  double                      EIR;  /* daily EIR */
};

#endif
