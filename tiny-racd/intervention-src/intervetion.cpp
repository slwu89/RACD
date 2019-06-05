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

intervention_manager::intervention_manager(house_vector* houses_) : houses(houses_) {};

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
