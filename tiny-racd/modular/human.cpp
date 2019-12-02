/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  April 2019
 #
 #  the human
*/

#include "human.hpp"





/* ################################################################################
#   track output
################################################################################ */

Rcpp::NumericVector track_transmission(human const* const human_ptr){
  return Rcpp::NumericVector::create(
    Rcpp::Named("epsilon") = human_ptr->epsilon,
    Rcpp::Named("lambda") = human_ptr->lambda,
    Rcpp::Named("phi") = human_ptr->phi,
    Rcpp::Named("prDetectAMic") = human_ptr->prDetectAMic,
    Rcpp::Named("prDetectAPCR") = human_ptr->prDetectAPCR,
    Rcpp::Named("prDetectUPCR") = human_ptr->prDetectUPCR,
    Rcpp::Named("c") = human_ptr->c
  );
};

Rcpp::NumericVector track_immunity(human const* const human_ptr){
  return Rcpp::NumericVector::create(
    Rcpp::Named("IB") = human_ptr->IB,
    Rcpp::Named("ID") = human_ptr->ID,
    Rcpp::Named("ICA") = human_ptr->ICA,
    Rcpp::Named("ICM") = human_ptr->ICM
  );
};

Rcpp::List human_2list(human const* const human_ptr){
  return Rcpp::List::create(
    Rcpp::Named("id") = human_ptr->id,
    Rcpp::Named("age") = human_ptr->age,
    Rcpp::Named("zeta") = human_ptr->zeta,
    Rcpp::Named("IB") = human_ptr->IB,
    Rcpp::Named("ID") = human_ptr->ID,
    Rcpp::Named("ICA") = human_ptr->ICA,
    Rcpp::Named("ICM") = human_ptr->ICM,
    Rcpp::Named("epsilon") = human_ptr->epsilon,
    Rcpp::Named("lambda") = human_ptr->lambda,
    Rcpp::Named("phi") = human_ptr->phi,
    Rcpp::Named("prDetectAMic") = human_ptr->prDetectAMic,
    Rcpp::Named("prDetectAPCR") = human_ptr->prDetectAPCR,
    Rcpp::Named("prDetectUPCR") = human_ptr->prDetectUPCR,
    Rcpp::Named("c") = human_ptr->c,
    Rcpp::Named("state") = human_ptr->state
  );
};
