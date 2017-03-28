// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// C_dense_dist
NumericMatrix C_dense_dist(std::string dist_type, const IntegerVector& XIX, const NumericMatrix& xmat, const IntegerVector& YIX, const NumericMatrix& ymat, bool pairwise, bool self);
RcppExport SEXP simdist_C_dense_dist(SEXP dist_typeSEXP, SEXP XIXSEXP, SEXP xmatSEXP, SEXP YIXSEXP, SEXP ymatSEXP, SEXP pairwiseSEXP, SEXP selfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type dist_type(dist_typeSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type XIX(XIXSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type xmat(xmatSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type YIX(YIXSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type ymat(ymatSEXP);
    Rcpp::traits::input_parameter< bool >::type pairwise(pairwiseSEXP);
    Rcpp::traits::input_parameter< bool >::type self(selfSEXP);
    rcpp_result_gen = Rcpp::wrap(C_dense_dist(dist_type, XIX, xmat, YIX, ymat, pairwise, self));
    return rcpp_result_gen;
END_RCPP
}
// C_sparse_dist
NumericMatrix C_sparse_dist(std::string dist_type, IntegerVector& XIX, IntegerVector& XP, NumericVector& XV, IntegerVector& YIX, IntegerVector& YP, NumericVector& YV, const size_t out_rows, const size_t out_cols, const bool pairwise, const bool self, const int xreorder, const int yreorder);
RcppExport SEXP simdist_C_sparse_dist(SEXP dist_typeSEXP, SEXP XIXSEXP, SEXP XPSEXP, SEXP XVSEXP, SEXP YIXSEXP, SEXP YPSEXP, SEXP YVSEXP, SEXP out_rowsSEXP, SEXP out_colsSEXP, SEXP pairwiseSEXP, SEXP selfSEXP, SEXP xreorderSEXP, SEXP yreorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type dist_type(dist_typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type XIX(XIXSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type XP(XPSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XV(XVSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type YIX(YIXSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type YP(YPSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type YV(YVSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_rows(out_rowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_cols(out_colsSEXP);
    Rcpp::traits::input_parameter< const bool >::type pairwise(pairwiseSEXP);
    Rcpp::traits::input_parameter< const bool >::type self(selfSEXP);
    Rcpp::traits::input_parameter< const int >::type xreorder(xreorderSEXP);
    Rcpp::traits::input_parameter< const int >::type yreorder(yreorderSEXP);
    rcpp_result_gen = Rcpp::wrap(C_sparse_dist(dist_type, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, pairwise, self, xreorder, yreorder));
    return rcpp_result_gen;
END_RCPP
}
// C_triplet_dist
NumericMatrix C_triplet_dist(std::string dist_type, IntegerVector& xpri, IntegerVector& xsec, NumericVector& xval, IntegerVector& ypri, IntegerVector& ysec, NumericVector& yval, const size_t out_rows, const size_t out_cols, const bool pairwise, const bool self, const int xreorder, const int yreorder);
RcppExport SEXP simdist_C_triplet_dist(SEXP dist_typeSEXP, SEXP xpriSEXP, SEXP xsecSEXP, SEXP xvalSEXP, SEXP ypriSEXP, SEXP ysecSEXP, SEXP yvalSEXP, SEXP out_rowsSEXP, SEXP out_colsSEXP, SEXP pairwiseSEXP, SEXP selfSEXP, SEXP xreorderSEXP, SEXP yreorderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type dist_type(dist_typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type xpri(xpriSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type xsec(xsecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type xval(xvalSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ypri(ypriSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ysec(ysecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yval(yvalSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_rows(out_rowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_cols(out_colsSEXP);
    Rcpp::traits::input_parameter< const bool >::type pairwise(pairwiseSEXP);
    Rcpp::traits::input_parameter< const bool >::type self(selfSEXP);
    Rcpp::traits::input_parameter< const int >::type xreorder(xreorderSEXP);
    Rcpp::traits::input_parameter< const int >::type yreorder(yreorderSEXP);
    rcpp_result_gen = Rcpp::wrap(C_triplet_dist(dist_type, xpri, xsec, xval, ypri, ysec, yval, out_rows, out_cols, pairwise, self, xreorder, yreorder));
    return rcpp_result_gen;
END_RCPP
}
// C_sparse_aggr_dist
NumericMatrix C_sparse_aggr_dist(const std::string& aggr_type, const std::string& dist_type, NumericMatrix& vecs, IntegerVector& XIX, IntegerVector& XP, NumericVector& XV, IntegerVector& YIX, IntegerVector& YP, NumericVector& YV, const size_t out_rows, const size_t out_cols, const bool pairwise, const bool self, const bool precompute);
RcppExport SEXP simdist_C_sparse_aggr_dist(SEXP aggr_typeSEXP, SEXP dist_typeSEXP, SEXP vecsSEXP, SEXP XIXSEXP, SEXP XPSEXP, SEXP XVSEXP, SEXP YIXSEXP, SEXP YPSEXP, SEXP YVSEXP, SEXP out_rowsSEXP, SEXP out_colsSEXP, SEXP pairwiseSEXP, SEXP selfSEXP, SEXP precomputeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type aggr_type(aggr_typeSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type dist_type(dist_typeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type vecs(vecsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type XIX(XIXSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type XP(XPSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type XV(XVSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type YIX(YIXSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type YP(YPSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type YV(YVSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_rows(out_rowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_cols(out_colsSEXP);
    Rcpp::traits::input_parameter< const bool >::type pairwise(pairwiseSEXP);
    Rcpp::traits::input_parameter< const bool >::type self(selfSEXP);
    Rcpp::traits::input_parameter< const bool >::type precompute(precomputeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_sparse_aggr_dist(aggr_type, dist_type, vecs, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, pairwise, self, precompute));
    return rcpp_result_gen;
END_RCPP
}
// C_triplet_aggr_dist
NumericMatrix C_triplet_aggr_dist(const std::string& aggr_type, const std::string& dist_type, NumericMatrix& vecs, IntegerVector& xpri, IntegerVector& xsec, NumericVector& xval, IntegerVector& ypri, IntegerVector& ysec, NumericVector& yval, const size_t out_rows, const size_t out_cols, const bool pairwise, const bool self, const bool precompute);
RcppExport SEXP simdist_C_triplet_aggr_dist(SEXP aggr_typeSEXP, SEXP dist_typeSEXP, SEXP vecsSEXP, SEXP xpriSEXP, SEXP xsecSEXP, SEXP xvalSEXP, SEXP ypriSEXP, SEXP ysecSEXP, SEXP yvalSEXP, SEXP out_rowsSEXP, SEXP out_colsSEXP, SEXP pairwiseSEXP, SEXP selfSEXP, SEXP precomputeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type aggr_type(aggr_typeSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type dist_type(dist_typeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type vecs(vecsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type xpri(xpriSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type xsec(xsecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type xval(xvalSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ypri(ypriSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ysec(ysecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type yval(yvalSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_rows(out_rowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type out_cols(out_colsSEXP);
    Rcpp::traits::input_parameter< const bool >::type pairwise(pairwiseSEXP);
    Rcpp::traits::input_parameter< const bool >::type self(selfSEXP);
    Rcpp::traits::input_parameter< const bool >::type precompute(precomputeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_triplet_aggr_dist(aggr_type, dist_type, vecs, xpri, xsec, xval, ypri, ysec, yval, out_rows, out_cols, pairwise, self, precompute));
    return rcpp_result_gen;
END_RCPP
}
