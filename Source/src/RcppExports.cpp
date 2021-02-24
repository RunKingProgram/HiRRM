// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// calynew
arma::mat calynew(arma::mat ynew, arma::mat vl);
RcppExport SEXP _HiRRM_calynew(SEXP ynewSEXP, SEXP vlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ynew(ynewSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vl(vlSEXP);
    rcpp_result_gen = Rcpp::wrap(calynew(ynew, vl));
    return rcpp_result_gen;
END_RCPP
}
// fast_mvlmmthreads
arma::mat fast_mvlmmthreads(int thread, arma::vec Y, arma::mat gnew, arma::mat vl, arma::vec sth);
RcppExport SEXP _HiRRM_fast_mvlmmthreads(SEXP threadSEXP, SEXP YSEXP, SEXP gnewSEXP, SEXP vlSEXP, SEXP sthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type thread(threadSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gnew(gnewSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type vl(vlSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sth(sthSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_mvlmmthreads(thread, Y, gnew, vl, sth));
    return rcpp_result_gen;
END_RCPP
}
// MatMult
SEXP MatMult(arma::mat A, arma::mat B);
RcppExport SEXP _HiRRM_MatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenMatMult
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _HiRRM_eigenMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenMapMatMult
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _HiRRM_eigenMapMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMapMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// CLegendre
Eigen::MatrixXd CLegendre(Eigen::VectorXd x, int s);
RcppExport SEXP _HiRRM_CLegendre(SEXP xSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(CLegendre(x, s));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HiRRM_calynew", (DL_FUNC) &_HiRRM_calynew, 2},
    {"_HiRRM_fast_mvlmmthreads", (DL_FUNC) &_HiRRM_fast_mvlmmthreads, 5},
    {"_HiRRM_MatMult", (DL_FUNC) &_HiRRM_MatMult, 2},
    {"_HiRRM_eigenMatMult", (DL_FUNC) &_HiRRM_eigenMatMult, 2},
    {"_HiRRM_eigenMapMatMult", (DL_FUNC) &_HiRRM_eigenMapMatMult, 2},
    {"_HiRRM_CLegendre", (DL_FUNC) &_HiRRM_CLegendre, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_HiRRM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
