#include "common.h"
#include "distances.h"
#include "workers_dist.h"
#include "workers_aggr.h"

using namespace simdist;

// [[Rcpp::export]]
NumericMatrix C_dense_dist(std::string dist_type,
                           const IntegerVector& XIX, const NumericMatrix& xmat,
                           const IntegerVector& YIX, const NumericMatrix& ymat,
                           bool pairwise = false, bool self = false) {
  if (dist_type == "COSINE")
    return denseDist<COSINE>(XIX, xmat, YIX, ymat, pairwise, self);
  else if (dist_type == "EUCLIDEAN")
    return denseDist<EUCLIDEAN>(XIX, xmat, YIX, ymat, pairwise, self);
  else
    throw std::invalid_argument("Invalid dist_type: " + dist_type);
}

/* // [[Rcpp::export]] */
/* NumericMatrix C_dense_range_dist(std::string dist_type, */
/*                                  const IntegerVector& xrange, const NumericMatrix& XV, */
/*                                  const IntegerVector& yrange, const NumericMatrix& YV, */
/*                                  bool self = false) { */
/*   if (dist_type == "COSINE") */
/*     return denseRangeDist<COSINE>(xrange, XV, yrange, YV, self); */
/*   else if (dist_type == "EUCLIDEAN") */
/*     return denseRangeDist<EUCLIDEAN>(xrange, XV, yrange, YV, self); */
/*   else */
/*     throw std::invalid_argument("Invalid dist_type: " + dist_type); */
/* } */

// [[Rcpp::export]]
NumericMatrix C_sparse_dist(std::string dist_type,
                            IntegerVector& XIX, IntegerVector& XP, NumericVector& XV,
                            IntegerVector& YIX, IntegerVector& YP, NumericVector& YV,
                            const size_t out_rows, const size_t out_cols,
                            const bool pairwise = false, const bool self = false,
                            const int xreorder = 0, const int yreorder = 0) {
  if (dist_type == "COSINE") 
    return sparseDist<COSINE>(XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, pairwise, self, xreorder, yreorder);
  else if (dist_type == "EUCLIDEAN")
    return sparseDist<EUCLIDEAN>(XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, pairwise, self, xreorder, yreorder);
  else
    throw std::invalid_argument("Invalid dist_type: " + dist_type);
}

// [[Rcpp::export]]
NumericMatrix C_triplet_dist(std::string dist_type,
                             IntegerVector& xpri, IntegerVector& xsec, NumericVector& xval,
                             IntegerVector& ypri, IntegerVector& ysec, NumericVector& yval,
                             const size_t out_rows, const size_t out_cols,
                             const bool pairwise = false, const bool self = false,
                             const int xreorder = 0, const int yreorder = 0) {
  if (dist_type == "COSINE") 
    return tripletDist<COSINE>(xpri, xsec, xval, ypri, ysec, yval, out_rows, out_cols, pairwise, self, xreorder, yreorder);
  else if (dist_type == "EUCLIDEAN")
    return tripletDist<EUCLIDEAN>(xpri, xsec, xval, ypri, ysec, yval, out_rows, out_cols, pairwise, self, xreorder, yreorder);
  else
    throw std::invalid_argument("Invalid dist_type: " + dist_type);
}


// [[Rcpp::export]]
NumericMatrix C_sparse_aggr_dist(const std::string& aggr_type, const std::string& dist_type,
                                 NumericMatrix& vecs,
                                 IntegerVector& XIX, IntegerVector& XP,
                                 NumericVector& XV, IntegerVector& YIX,
                                 IntegerVector& YP, NumericVector& YV,
                                 const size_t out_rows, const size_t out_cols,
                                 const bool self = false, const bool precompute = true) {
  if (aggr_type == "CENTROID") {
    if (dist_type == "COSINE") {
      return sparseAggrDist<CENTROID, COSINE>(vecs, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, self, precompute);
    } else if (dist_type == "EUCLIDEAN") {
      return sparseAggrDist<CENTROID, EUCLIDEAN>(vecs, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, self, precompute);
    } else {
      throw std::invalid_argument("Invalid dist_type: " + dist_type);
    }
  } else if (aggr_type == "SEM_MIN_MAX") {
    if (dist_type == "COSINE") {
      return sparseAggrDist<SEM_MIN_MAX, COSINE>(vecs, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, self, precompute);
    } else if (dist_type == "EUCLIDEAN") {
      return sparseAggrDist<SEM_MIN_MAX, EUCLIDEAN>(vecs, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, self, precompute);
    } else {
      throw std::invalid_argument("Invalid dist_type: " + dist_type);
    }
  } else if (aggr_type == "SEM_MIN_SUM") {
    if (dist_type == "COSINE") {
      return sparseAggrDist<SEM_MIN_SUM, COSINE>(vecs, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, self, precompute);
    } else if (dist_type == "EUCLIDEAN") {
      return sparseAggrDist<SEM_MIN_SUM, EUCLIDEAN>(vecs, XIX, XP, XV, YIX, YP, YV, out_rows, out_cols, self, precompute);
    } else {
      throw std::invalid_argument("Invalid dist_type: " + dist_type);
    }
  } else {
    throw std::invalid_argument("Invalid aggr_type: " + aggr_type);
  }
}

// [[Rcpp::export]]
NumericMatrix C_triplet_aggr_dist(const std::string& aggr_type, const std::string& dist_type,
                                  NumericMatrix& vecs,
                                  IntegerVector& xpri, IntegerVector& xsec, NumericVector& xval,
                                  IntegerVector& ypri, IntegerVector& ysec, NumericVector& yval,
                                  const size_t out_rows, const size_t out_cols,
                                  const bool self = false, const bool precompute = true) {
  auto X = TripletMat<IntegerVector, NumericVector>(xpri, xsec, xval, out_rows).toSparseMat();
  auto Y = TripletMat<IntegerVector, NumericVector>(ypri, ysec, yval, out_cols).toSparseMat();
  return C_sparse_aggr_dist(aggr_type, dist_type, vecs,
                            X.i, X.p, X.v, Y.i, Y.p, Y.v,
                            out_rows, out_cols,
                            self, precompute);
}
