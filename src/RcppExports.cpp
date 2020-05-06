// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// d_snb
NumericVector d_snb(NumericVector& x, Nullable<NumericVector> size_param, Nullable<NumericVector> prob_param, Nullable<NumericVector> mu_param, const int& infinite, const bool& log_v);
RcppExport SEXP _stochprofML_d_snb(SEXP xSEXP, SEXP size_paramSEXP, SEXP prob_paramSEXP, SEXP mu_paramSEXP, SEXP infiniteSEXP, SEXP log_vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type size_param(size_paramSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type prob_param(prob_paramSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type mu_param(mu_paramSEXP);
    Rcpp::traits::input_parameter< const int& >::type infinite(infiniteSEXP);
    Rcpp::traits::input_parameter< const bool& >::type log_v(log_vSEXP);
    rcpp_result_gen = Rcpp::wrap(d_snb(x, size_param, prob_param, mu_param, infinite, log_v));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_stochprofML_d_snb", (DL_FUNC) &_stochprofML_d_snb, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_stochprofML(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
