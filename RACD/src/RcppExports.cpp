// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// test_parameters
void test_parameters(const Rcpp::NumericVector& theta);
RcppExport SEXP _RACD_test_parameters(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta(thetaSEXP);
    test_parameters(theta);
    return R_NilValue;
END_RCPP
}
// test_human_parameters
void test_human_parameters(const Rcpp::NumericVector& theta);
RcppExport SEXP _RACD_test_human_parameters(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta(thetaSEXP);
    test_human_parameters(theta);
    return R_NilValue;
END_RCPP
}
// test_prng
void test_prng(const uint_least32_t& seed);
RcppExport SEXP _RACD_test_prng(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const uint_least32_t& >::type seed(seedSEXP);
    test_prng(seed);
    return R_NilValue;
END_RCPP
}
// test_house
void test_house(const Rcpp::NumericVector& theta, const uint_least32_t& seed);
RcppExport SEXP _RACD_test_house(SEXP thetaSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const uint_least32_t& >::type seed(seedSEXP);
    test_house(theta, seed);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RACD_test_parameters", (DL_FUNC) &_RACD_test_parameters, 1},
    {"_RACD_test_human_parameters", (DL_FUNC) &_RACD_test_human_parameters, 1},
    {"_RACD_test_prng", (DL_FUNC) &_RACD_test_prng, 1},
    {"_RACD_test_house", (DL_FUNC) &_RACD_test_house, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RACD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}