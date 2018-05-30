#include "RACD-Mosquito.hpp"

/* ################################################################################
 * abstract base class
################################################################################ */

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


/* ################################################################################
 * default class (seasonally forced EIR)
################################################################################ */

/* constructor */
mosquito_habitat_eir::mosquito_habitat_eir(const double EIR_mean_, const double offset_, const int habitatID_, village* village_ptr_) :
  mosquito_habitat_base(habitatID_,village_ptr_), EIR_mean(EIR_mean_), offset(offset_) {

    #ifdef DEBUG_RACD
    std::cout << "mosquito_habitat_eir " << habitatID << " being born at " << this << std::endl;
    #endif

};

/* destructor */
mosquito_habitat_eir::~mosquito_habitat_eir(){
  #ifdef DEBUG_RACD
  std::cout << "mosquito_habitat_eir " << habitatID << " being killed at " << this << std::endl;
  #endif
};

/* daily simulation */
void mosquito_habitat_eir::one_day(const int tNow){
  EIR = EIR_mean * (1+sin(2*M_PI*((double)tNow-offset)/365.25));
};
