#include "RACD-MosquitoHabitat.hpp"

/* abstract base class */

/* constructor for base class */
mosquito_habitat_base::mosquito_habitat_base(const int habitatID_, village* village_ptr_) :
  habitatID(habitatID_), village_ptr(village_ptr_) {
  #ifdef DEBUG_RACD
  std::cout << "mosquito_habitat_base " << habitatID << " being born at " << this << std::endl;
  #endif
};

/* destructor */
mosquito_habitat_base::~mosquito_habitat_base(){
  #ifdef DEBUG_RACD
  std::cout << "mosquito_habitat_base " << habitatID << " being killed at " << this << std::endl;
  #endif
}
