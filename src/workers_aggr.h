#ifndef DISTANCE_WORKERS_AGGR_HPP__
#define DISTANCE_WORKERS_AGGR_HPP__

#include "common.h"
#include "workers_dist.h"       // for denseDist
#include "workers_misc.h" // centroid and normalization

namespace simdist {

 /// BASE WORKER

// Worker that takes two sparse matrices and computes output distance matrix
// element by element
template<DistType DT>
struct SparseAggrBaseWorker : public RcppParallel::Worker {
  
  const Dist<DT> dist {};
  const bool self;

  IntegerVector XIX0;
  RVector<int> XIX;
  const RVector<int> XP;
  const RVector<double> XV;

  IntegerVector YIX0;
  RVector<int> YIX;
  const RVector<int> YP;
  const RVector<double> YV;

  NumericMatrix vecs0;
  RMatrix<double> vecs;
  RMatrix<double> outmat;

  explicit SparseAggrBaseWorker(const NumericMatrix& vecs,
                                IntegerVector& XIX,
                                const IntegerVector& XP,
                                const NumericVector& XV,
                                IntegerVector& YIX,
                                const IntegerVector& YP,
                                const NumericVector& YV,
                                NumericMatrix& outmat,
                                bool self = false)
    : vecs0(vecs), vecs(vecs), 
      XIX(XIX), XIX0(XIX), XP(XP), XV(XV),
      YIX(YIX), YIX0(YIX), YP(YP), YV(YV), 
      outmat(outmat), self(self) { }

  void operator()(std::size_t xjbeg, std::size_t xjend) {

    for (size_t xj = xjbeg; xj < xjend; xj++) {

      size_t yjend = this->self ? xj : this->outmat.ncol();

      for (std::size_t yj = 0; yj < yjend; yj++) {
        int xbeg = this->XP[xj], xend = this->XP[xj + 1];
        int ybeg = this->YP[yj], yend = this->YP[yj + 1];


        if (xbeg == xend || ybeg == yend)  {

          this->outmat(xj, yj) = std::numeric_limits<double>::quiet_NaN();

        } else {

          this->outmat(xj, yj) = computeOne(xj, xbeg, xend, yj, ybeg, yend);        
        }

        if (self)
          this->outmat(yj, xj) = this->outmat(xj, yj);
      }
    }
  }

  virtual void precompute() {};
  virtual double computeOne(const size_t xix, const int xbeg, const int xend,
                            const size_t yix, const int ybeg, const int yend) = 0;
     
};



 /// COMMON 3D DIST TEMPLATE

template<AggrType AT, DistType DT>
struct SparseAggrWorker : public SparseAggrBaseWorker<DT> {};



 /// CENTROID DISTANCE

template<DistType DT>
struct SparseAggrWorker<CENTROID, DT> : public SparseAggrBaseWorker<DT>
{
  bool precomputed = false;
  NumericMatrix Xvecs, Yvecs;
  
  // bring constructor into this scope
  using SparseAggrBaseWorker<DT>::SparseAggrBaseWorker;

  void precompute()
  {
    precomputed = true;
    Xvecs = NumericMatrix(this->vecs.nrow(), this->outmat.nrow());
    CentroidWorker xworker(this->vecs, this->XIX, this->XP, this->XV, Xvecs);
    parallelFor(0, this->outmat.nrow(), xworker);
    
    Yvecs = NumericMatrix(this->vecs.nrow(), this->outmat.ncol());
    CentroidWorker yworker(this->vecs, this->YIX, this->YP, this->YV, Yvecs);
    parallelFor(0, this->outmat.ncol(), yworker);

    // hardcoding L2 normalization for cosine distance
    if (DT == COSINE) {
      DenseNormWorker<L2> xnormer(Xvecs);
      parallelFor(0, Xvecs.ncol(), xnormer, 10);
      DenseNormWorker<L2> ynormer(Yvecs);
      parallelFor(0, Yvecs.ncol(), ynormer, 10);
    }    
  }
  
  double computeOne(const size_t xix, const int xbeg, const int xend,
                    const size_t yix, const int ybeg, const int yend)
  {
    if (precomputed) {
      return this->dist(Xvecs.column(xix).begin(), Xvecs.column(xix).end(), Yvecs.column(yix).begin());
    } else {
      throw std::runtime_error("precomputed=FALSE is not yet implemented for Centroid Aggregation Distance");
    }
  }
};



 /// RWMD (RELAXED WORD MOVER DISTANCE)

namespace {

// IX is modified by side effect, return new range
inline IntegerVector compact_range(RVector<int>& IX, const size_t nvecs) {

  std::unordered_map<int, int> ixmap;
  std::vector<int> range;

  for (int& ix : IX) {
    if (ix == NA_INTEGER)
      throw std::out_of_range("NA indexing not allowed");
    auto ixnew_it = ixmap.find(ix);
    if (ixnew_it == ixmap.end()) {
      if (ix >= nvecs)
        throw std::out_of_range("Index" + std::to_string(ix) +
                                " larger than number of word vectors");
      int ixnew = ixmap.size();
      ixmap[ix] = ixnew;
      range.push_back(ix);
      ix = ixnew;
    } else {
      ix = ixnew_it->second;
    }
  }
    
  return IntegerVector(range.begin(), range.end());
}

}


template<DistType DT>
struct SparseAggrMinWorker : public SparseAggrBaseWorker<DT> {

  bool precomputed = false;
  NumericMatrix precomputed_dist;
  
  using SparseAggrBaseWorker<DT>::SparseAggrBaseWorker;

  virtual double aggr_fun(double xmin, double ymin) = 0;

  void precompute() {

    precomputed = true;
    
    if (this->self) {
      this->XIX0 = clone(this->XIX0);
      this->XIX = RVector<int>(this->XIX0);
      IntegerVector range = compact_range(this->XIX, this->vecs.ncol());
      precomputed_dist = denseRangeDist<DT>(range, range, this->vecs0, true);
    } else {
      this->XIX0 = clone(this->XIX0);
      this->XIX = RVector<int>(this->XIX0);
      this->YIX0 = clone(this->YIX0);
      this->YIX = RVector<int>(this->YIX0);
      IntegerVector xrange = compact_range(this->XIX, this->vecs.ncol());
      IntegerVector yrange = compact_range(this->YIX, this->vecs.ncol());
      precomputed_dist = denseRangeDist<DT>(xrange, yrange, this->vecs0, false);
    }
  }
  
  double computeOne(const size_t xj, const int xbeg, const int xend,
                    const size_t yj, const int ybeg, const int yend) {

    int xlen = xend - xbeg, ylen = yend - ybeg;
    std::vector<double> xmins(xlen, POS_INF), ymins(ylen, POS_INF);
    
    for (int r = 0; r < xlen; r++) {
      int xix = this->XIX[xbeg + r];
      for (int c = 0; c < ylen; c++) {
        int yix = this->YIX[ybeg + c];
        double dst;
        if (precomputed) {
          dst = precomputed_dist(xix, yix);
        } else {
          auto xcol = this->vecs.column(xix);
          auto ycol = this->vecs.column(yix);
          dst = this->dist(xcol.begin(), xcol.end(), ycol.begin());
        }
        xmins[r] = std::min(dst, xmins[r]);
        ymins[c] = std::min(dst, ymins[c]);
      }
    }

    return aggr_fun(std::inner_product(xmins.begin(), xmins.end(), this->XV.begin() + xbeg, 0.0),
                    std::inner_product(ymins.begin(), ymins.end(), this->YV.begin() + ybeg, 0.0));
  }
  
};

template<DistType DT>
struct SparseAggrWorker<SEM_MIN_MAX, DT> : public SparseAggrMinWorker<DT>
{
  using SparseAggrMinWorker<DT>::SparseAggrMinWorker;
  inline double aggr_fun(double xmin, double ymin) {
    return std::max(xmin, ymin);
  }
};
  
template<DistType DT>
struct SparseAggrWorker<SEM_MIN_SUM, DT> : public SparseAggrMinWorker<DT>
{
  using SparseAggrMinWorker<DT>::SparseAggrMinWorker;
  inline double aggr_fun(double xmin, double ymin) {
    return xmin + ymin;
  }
};



// ENTRY POINT

template<AggrType AT, DistType DT>
NumericMatrix sparseAggrDist(NumericMatrix& vecs,
                             IntegerVector& XIX, IntegerVector& XP, NumericVector& XV,
                             IntegerVector& YIX, IntegerVector& YP, NumericVector& YV,
                             const size_t out_rows, const size_t out_cols,
                             const bool self = false, const bool precompute = true)
{
  NumericMatrix outmat(out_rows, out_cols);
  SparseAggrWorker<AT, DT> worker(vecs, XIX, XP, XV, YIX, YP, YV, outmat, false);
  if (precompute) worker.precompute();
  parallelFor(0, out_rows, worker);
  return outmat;
}


} // namespace simdist

#endif
