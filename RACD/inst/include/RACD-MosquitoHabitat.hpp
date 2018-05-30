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
#include <memory>

/* alias and forward declarations */
class village;

/*
 * abstract base class
 *  basically, the mosquito habitats need to simulate on a one day time step
 *  and send out bites to houses. They have a pointer to the containing village
 *  that allows them to do this.
 *
*/
class mosquito_habitat_base {

public:
  /* define the interface */

  /* constructor */
  mosquito_habitat_base(const int habitatID_, village* village_ptr_);

  /* pure virtual destructor */
  virtual ~mosquito_habitat_base() = 0;

  /* required implementations */
  virtual void one_day() = 0;
  virtual void distribute_bites() = 0;

private:
  /* default data members */
  int                         habitatID;
  village*                    village_ptr; /* raw pointer ok because house lifespan > human lifespan in memory */
};

#endif
