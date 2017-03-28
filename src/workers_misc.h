// -*- mode: c++ -*-

#ifndef DISTANCE_WORKERS_MISC_HPP__
#define DISTANCE_WORKERS_MISC_HPP__

#include "norms.h"
#include "common.h"

namespace simdist {

 /// DENSE COLUMN WORKERS

template<NormType NT>
struct DenseNormWorker : public RcppParallel::Worker {

  Norm<NT> norm {};
  RMatrix<double> outmat;

  explicit DenseNormWorker(NumericMatrix& outmat) : outmat(outmat) { }
  explicit DenseNormWorker(RMatrix<double>& outmat) : outmat(outmat) { }

  void operator()(std::size_t xbeg, std::size_t xend) {
    for (size_t xj = xbeg; xj < xend; xj++) {
      norm(this->outmat.column(xj).begin(), this->outmat.column(xj).end());
    }
  }
};

struct DenseColumnWorker : public RcppParallel::Worker {
  
  RMatrix<double> vecs;
  RMatrix<double> inmat;
  RMatrix<double> outmat;

  explicit DenseColumnWorker(NumericMatrix& vecs,
                             NumericMatrix& inmat,
                             NumericMatrix& outmat)
    : vecs(vecs), inmat(inmat), outmat(outmat) { }

  explicit DenseColumnWorker(RMatrix<double>& vecs,
                             NumericMatrix& inmat, 
                             NumericMatrix& outmat)
    : vecs(vecs), inmat(inmat), outmat(outmat) { }

  void operator()(std::size_t xbeg, std::size_t xend) {
    for (size_t xj = xbeg; xj < xend; xj++) {
      computeColumn(xj);
    }
  }

  virtual void computeColumn(const size_t xix) = 0;
};


 /// SPARSE COLUMN WORKERS

// Worker that takes one sparse matrix as input and computes output matrix
// column by column
struct SparseColumnWorker : public RcppParallel::Worker {
  
  RVector<int> IX;
  const RVector<int> P;
  const RVector<double> V;
  RMatrix<double> vecs;
  RMatrix<double> outmat;

  explicit SparseColumnWorker(const NumericMatrix& vecs,
                              IntegerVector& IX,
                              const IntegerVector& P,
                              const NumericVector& V,
                              NumericMatrix& outmat)
    : vecs(vecs), IX(IX), P(P), V(V), outmat(outmat) { }

  explicit SparseColumnWorker(const RMatrix<double>& vecs,
                              RVector<int>& IX,
                              const RVector<int>& P,
                              const RVector<double>& V,
                              NumericMatrix& outmat)
    : vecs(vecs), IX(IX), P(P), V(V), outmat(outmat) { }

  void operator()(std::size_t xbeg, std::size_t xend) {
    for (size_t xj = xbeg; xj < xend; xj++) {
      computeColumn(xj, P[xj], P[xj + 1]);
    }
  }

  virtual void computeColumn(const size_t xix, const int xbeg, const int xend) = 0;
};

// Compute mean vector for each xix
struct CentroidWorker : public SparseColumnWorker
{
  using SparseColumnWorker::SparseColumnWorker;
  
  void computeColumn(const size_t xix, const int xbeg, const int xend)
  {
    size_t nrows = this->vecs.nrow();
    vector<double> tt(nrows);

    for (int xi = xbeg; xi < xend; xi++) {
      auto vec = this->vecs.column(this->IX[xi]);
      double val = V[xi];
      for (size_t vi = 0;  vi < nrows; vi++) {
        tt[vi] +=  val * vec[vi];
      }
    }

    size_t nvecs = xend - xbeg;
    for (size_t i = 0; i < nrows; i++) {
      this->outmat(i, xix) = tt[i]/nvecs;
    }
  }  
};

} // namespace simdist

#endif
