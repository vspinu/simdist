#ifndef DISTANCE_WORKERS_HPP__
#define DISTANCE_WORKERS_HPP__

#include "distances.h"
#include "sparse_mat.h"

namespace simdist {


/// DENSE WORKERS

template <DistType DT>
struct DenseMatWorker : public RcppParallel::Worker {

  const Dist<DT> dist {};
  const bool self;
  const RMatrix<double> xmat; 
  const RMatrix<double> ymat; 
  RMatrix<double> outmat;

  DenseMatWorker(const NumericMatrix& xmat,
                 const NumericMatrix& ymat,
                 NumericMatrix& outmat,
                 bool self = false)
    : xmat(xmat), ymat(ymat), outmat(outmat), self(self) {}

  void operator()(std::size_t xbegin, std::size_t xend) {
    for (std::size_t xj = xbegin; xj < xend; xj++) {
      std::size_t yend = self ? xj : ymat.ncol();
      const auto& xcol = xmat.column(xj);
      for (std::size_t yj = 0; yj < yend; yj++) {
        const auto& ycol = ymat.column(yj);
        double dst = dist(xcol.begin(), xcol.end(), ycol.begin());
        outmat(xj, yj) = dst;
        if (self)
          outmat(yj, xj) = dst;
      }
    }
  }
};

template <DistType DT>
NumericMatrix denseDist(const NumericMatrix& xmat, const NumericMatrix& ymat) {
  if (xmat.nrow() != ymat.nrow())
    throw std::invalid_argument("xmat and ymat must have same number of rows");
  NumericMatrix outmat(xmat.ncol(), ymat.ncol());
  DenseMatWorker<DT> worker(xmat, ymat, outmat);
  parallelFor(0, xmat.ncol(), worker);
  return outmat;
}

// YIX and XIX are assumed to be sorted in ascending order
template <DistType DT>
struct DenseWorker : public RcppParallel::Worker {

  const Dist<DT> dist {};
  const bool self;
  const RVector<int> XIX;
  const RMatrix<double> xmat; 
  const RVector<int> YIX;
  const RMatrix<double> ymat; 
  RMatrix<double> outmat;

  DenseWorker(const IntegerVector& XIX,
              const NumericMatrix& xmat,
              const IntegerVector& YIX,
              const NumericMatrix& ymat,
              NumericMatrix& outmat,
              bool self = false)
    : XIX(XIX), xmat(xmat), YIX(YIX), ymat(ymat), outmat(outmat), self(self) {}

  void operator()(std::size_t xbegin, std::size_t xend) {
    for (std::size_t xj = xbegin; xj < xend; xj++) {
      std::size_t yend = self ? xj : ymat.ncol();
      const auto& xcol = xmat.column(xj);
      for (std::size_t yj = 0; yj < yend; yj++) {
        const auto& ycol = ymat.column(yj);
        double dst = dist(XIX.begin(), XIX.end(), xcol.begin(),
                          YIX.begin(), YIX.end(), ycol.begin());
        outmat(xj, yj) = dst;
        if (self)
          outmat(yj, xj) = dst;
      }
    }
  }
  
};

template <DistType DT>
NumericMatrix denseDist(const IntegerVector& XIX,
                        const NumericMatrix& xmat,
                        const IntegerVector& YIX,
                        const NumericMatrix& ymat,
                        const bool self = false) {
  NumericMatrix outmat(xmat.ncol(), ymat.ncol());
  DenseWorker<DT> worker(XIX, xmat, YIX, ymat, outmat, self);
  parallelFor(0, xmat.ncol(), worker);
  return outmat;
}

/// xrange and yrange are indexes in primary dimension
template <DistType DT>
struct DenseRangeWorker : public RcppParallel::Worker {

  const bool self;
  const RMatrix<double> xvecs;
  const RMatrix<double> yvecs;
  const RVector<int> xrange;
  const RVector<int> yrange;
  RMatrix<double> outmat;
  const Dist<DT> dist {};

  template<class VT>
  DenseRangeWorker(const VT& xrange,
                   const NumericMatrix& xvecs,
                   const VT& yrange,
                   const NumericMatrix& yvecs,
                   NumericMatrix& outmat,
                   bool self = false)
    : xrange(xrange), xvecs(xvecs), yrange(yrange), yvecs(yvecs), outmat(outmat), self(self) { }
  
  void operator()(std::size_t xbegin, std::size_t xend) {
    for (std::size_t xj = xbegin; xj < xend; xj++) {
      std::size_t yend = self ? xj : yrange.size();
      const auto& xcol = xvecs.column(xrange[xj]);
      for (std::size_t yj = 0; yj < yend; yj++) {
        const auto& ycol = yvecs.column(yrange[yj]);
        double dst = dist(xcol.begin(), xcol.end(), ycol.begin());
        outmat(xj, yj) = dst;
        if (self)
          outmat(yj, xj) = dst;
      }
    }
  }
  
};


template <DistType DT, class VT>
NumericMatrix denseRangeDist(const VT& xrange,
                             const VT& yrange,
                             const NumericMatrix& vecs,
                             const bool self = false) {
  NumericMatrix outmat(xrange.size(), yrange.size());
  DenseRangeWorker<DT> worker(xrange, vecs, yrange, vecs, outmat, self);
  parallelFor(0, xrange.size(), worker);
  return outmat;
}


template <DistType DT, class VT>
NumericMatrix denseRangeDist(const VT& xrange,
                             const NumericMatrix& xvecs,
                             const VT& yrange,
                             const NumericMatrix& yvecs,
                             const bool self = false) {
  NumericMatrix outmat(xrange.size(), yrange.size());
  DenseRangeWorker<DT> worker(xrange, xvecs, yrange, yvecs, outmat, self);
  parallelFor(0, xrange.size(), worker);
  return outmat;
}



/// SPARSE WORKERS

template<DistType DT>
struct SparseMatWorker : public RcppParallel::Worker {
  
  const Dist<DT> dist {};
  const bool self;

  const RVector<int> XIX;
  const RVector<int> XP;
  const RVector<double> XV;

  const RVector<int> YIX;
  const RVector<int> YP;
  const RVector<double> YV;

  RMatrix<double> outmat;
  bool precompute;
  NumericMatrix precomputed_mat;

  explicit SparseMatWorker(const IntegerVector& XIX, const IntegerVector& XP, const NumericVector& XV,
                           const IntegerVector& YIX, const IntegerVector& YP, const NumericVector& YV,
                           NumericMatrix outmat, const bool self = false)
    : XIX(XIX), XP(XP), XV(XV), YIX(YIX), YP(YP), YV(YV), 
      outmat(outmat), self(self) {}

  void operator()(std::size_t xbegin, std::size_t xend) {
    for (int xj = xbegin; xj < xend; xj++) {
      size_t yend = self ? xj : outmat.ncol();
      for (int yj = 0; yj < yend; yj++) {
        /* printf("xj:%d yj:%d  XP[xj]:%d  YP[yj]:%d XP[xj]:%d  YP[yj]:%d\n", */
        /*        xj, yj, XP[xj], YP[yj], XP[xj + 1], YP[yj + 1]); */
        double out = dist(XIX.begin() + XP[xj], XIX.begin() + XP[xj + 1], XV.begin() + XP[xj],
                          YIX.begin() + YP[yj], YIX.begin() + YP[yj + 1], YV.begin() + YP[yj]);
        outmat(xj, yj) = out;
        if (self) {
          outmat(yj, xj) = out;
        }
      }
    }
  }
  
};

template<DistType DT>
NumericMatrix sparseDist(IntegerVector& XIX, IntegerVector& XP, NumericVector& XV,
                         IntegerVector& YIX, IntegerVector& YP, NumericVector& YV,
                         const size_t out_rows, const size_t out_cols,
                         const bool self = false,
                         const int xreorder = 0, const int yreorder = 0) {
  if (xreorder == 2) {
    XIX = clone(XIX);
    XV = clone(XV);
  }
  if (yreorder == 2) {
    YIX = clone(YIX);
    YV = clone(YV);
  }
  if (xreorder) reorder(XIX, XP, XV);

  if (yreorder) {
    reorder(YIX, YP, YV);
  }

  NumericMatrix outmat(out_rows, out_cols);
  SparseMatWorker<DT> worker(XIX, XP, XV, YIX, YP, YV, outmat, self);
  parallelFor(0, out_rows, worker, 100);
  return outmat;
}

template<DistType DT>
NumericMatrix tripletDist(IntegerVector& xpri, IntegerVector& xsec, NumericVector& xval,
                          IntegerVector& ypri, IntegerVector& ysec, NumericVector& yval,
                          const size_t out_rows, const size_t out_cols, 
                          const bool self = false, int xreorder = 0, int yreorder = 0) {
  // always pass 1 as we are creating temporaries anyhow
  xreorder = xreorder ? 1 : 0;
  yreorder = yreorder ? 1 : 0;
  auto X = simdist::TripletMat<IntegerVector, NumericVector>(xpri, xsec, xval, out_rows).toSparseMat();
  auto Y = simdist::TripletMat<IntegerVector, NumericVector>(ypri, ysec, yval, out_cols).toSparseMat();
  return sparseDist<DT>(X.i, X.p, X.v, Y.i, Y.p, Y.v, out_rows, out_cols, self, xreorder, yreorder);
}

} // simdist namespace



#endif
