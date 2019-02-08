// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// C_calc_D
arma::vec C_calc_D(const double D, const int J, arma::rowvec pr0, Rcpp::List tpms);
RcppExport SEXP _openpopscr_C_calc_D(SEXP DSEXP, SEXP JSEXP, SEXP pr0SEXP, SEXP tpmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int >::type J(JSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type pr0(pr0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tpms(tpmsSEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_D(D, J, pr0, tpms));
    return rcpp_result_gen;
END_RCPP
}
// C_calc_llk
double C_calc_llk(const int n, const int J, const int M, const arma::mat pr0, const Rcpp::List pr_capture, const Rcpp::List tpms, const int num_cores, const int num_states, const arma::vec entry);
RcppExport SEXP _openpopscr_C_calc_llk(SEXP nSEXP, SEXP JSEXP, SEXP MSEXP, SEXP pr0SEXP, SEXP pr_captureSEXP, SEXP tpmsSEXP, SEXP num_coresSEXP, SEXP num_statesSEXP, SEXP entrySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type pr0(pr0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type pr_capture(pr_captureSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type tpms(tpmsSEXP);
    Rcpp::traits::input_parameter< const int >::type num_cores(num_coresSEXP);
    Rcpp::traits::input_parameter< const int >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type entry(entrySEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_llk(n, J, M, pr0, pr_capture, tpms, num_cores, num_states, entry));
    return rcpp_result_gen;
END_RCPP
}
// C_calc_pdet
double C_calc_pdet(const int J, arma::mat pr0, Rcpp::List pr_captures, Rcpp::List tpms, const int num_states);
RcppExport SEXP _openpopscr_C_calc_pdet(SEXP JSEXP, SEXP pr0SEXP, SEXP pr_capturesSEXP, SEXP tpmsSEXP, SEXP num_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type J(JSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pr0(pr0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pr_captures(pr_capturesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tpms(tpmsSEXP);
    Rcpp::traits::input_parameter< const int >::type num_states(num_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_pdet(J, pr0, pr_captures, tpms, num_states));
    return rcpp_result_gen;
END_RCPP
}
// C_calc_move_llk
double C_calc_move_llk(const int n, const int J, const arma::mat pr0, const Rcpp::List pr_capture, const Rcpp::List tpms, const arma::vec num_cells, const arma::vec inside, const double dx, const arma::vec dt, const arma::vec sd, const int num_cores, const int num_states, const arma::vec entry);
RcppExport SEXP _openpopscr_C_calc_move_llk(SEXP nSEXP, SEXP JSEXP, SEXP pr0SEXP, SEXP pr_captureSEXP, SEXP tpmsSEXP, SEXP num_cellsSEXP, SEXP insideSEXP, SEXP dxSEXP, SEXP dtSEXP, SEXP sdSEXP, SEXP num_coresSEXP, SEXP num_statesSEXP, SEXP entrySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type pr0(pr0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type pr_capture(pr_captureSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type tpms(tpmsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type num_cells(num_cellsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type inside(insideSEXP);
    Rcpp::traits::input_parameter< const double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< const int >::type num_cores(num_coresSEXP);
    Rcpp::traits::input_parameter< const int >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type entry(entrySEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_move_llk(n, J, pr0, pr_capture, tpms, num_cells, inside, dx, dt, sd, num_cores, num_states, entry));
    return rcpp_result_gen;
END_RCPP
}
// C_calc_move_pdet
double C_calc_move_pdet(const int J, arma::mat pr0, Rcpp::List pr_captures, Rcpp::List tpms, const arma::vec num_cells, const arma::vec inside, const double dx, const arma::vec dt, const arma::vec sd, const int num_states);
RcppExport SEXP _openpopscr_C_calc_move_pdet(SEXP JSEXP, SEXP pr0SEXP, SEXP pr_capturesSEXP, SEXP tpmsSEXP, SEXP num_cellsSEXP, SEXP insideSEXP, SEXP dxSEXP, SEXP dtSEXP, SEXP sdSEXP, SEXP num_statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type J(JSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pr0(pr0SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type pr_captures(pr_capturesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type tpms(tpmsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type num_cells(num_cellsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type inside(insideSEXP);
    Rcpp::traits::input_parameter< const double >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< const int >::type num_states(num_statesSEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_move_pdet(J, pr0, pr_captures, tpms, num_cells, inside, dx, dt, sd, num_states));
    return rcpp_result_gen;
END_RCPP
}
// C_calc_pr_capture
arma::field<arma::cube> C_calc_pr_capture(const int n, const int J, const int K, const int M, const arma::cube& capthist, const arma::cube& enc0, const arma::mat usage, const int num_cores, const int num_states, const int detector_type, const int n_prim, const arma::vec S, const arma::vec entry);
RcppExport SEXP _openpopscr_C_calc_pr_capture(SEXP nSEXP, SEXP JSEXP, SEXP KSEXP, SEXP MSEXP, SEXP capthistSEXP, SEXP enc0SEXP, SEXP usageSEXP, SEXP num_coresSEXP, SEXP num_statesSEXP, SEXP detector_typeSEXP, SEXP n_primSEXP, SEXP SSEXP, SEXP entrySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type J(JSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type capthist(capthistSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type enc0(enc0SEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type usage(usageSEXP);
    Rcpp::traits::input_parameter< const int >::type num_cores(num_coresSEXP);
    Rcpp::traits::input_parameter< const int >::type num_states(num_statesSEXP);
    Rcpp::traits::input_parameter< const int >::type detector_type(detector_typeSEXP);
    Rcpp::traits::input_parameter< const int >::type n_prim(n_primSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type S(SSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type entry(entrySEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_pr_capture(n, J, K, M, capthist, enc0, usage, num_cores, num_states, detector_type, n_prim, S, entry));
    return rcpp_result_gen;
END_RCPP
}
