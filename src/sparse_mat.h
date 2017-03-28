// -*- mode: c++ -*-

#ifndef SIMDIST_SPARSE_HPP__
#define SIMDIST_SPARSE_HPP__

namespace simdist {

template<class ixT, class valT>
struct TripletMat;

template<class ixT, class valT>
struct SparseMat;
  
template<class ixT, class valT>
struct TripletMat {
  ixT pix;
  ixT six;
  valT val;
  int psize;

  TripletMat(const ixT& pix, const ixT& six, const valT& val, const int& psize)
    : pix(pix), six(six), val(val), psize(psize)
  {
    if (pix.size() != six.size() || pix.size() != val.size()) {
      throw std::invalid_argument("pix, six and val must have same size");
    }
  }

  SparseMat<ixT, valT> toSparseMat() {

    ixT P(psize + 1);
    ixT IX(val.size());
    valT VAL(val.size());
    
    vector<vector<int>> xvv(psize);
    for (int i = 0; i <  pix.size(); i++) {
      int ix = pix[i];
      if (ix >= psize)
        throw std::out_of_range("primary indexes are not < psize");
      P[ix + 1]++;
      xvv[ix].push_back(i);
    }
    for (int i = 1; i < P.size(); i++)
      P[i] += P[i-1];

    int k = 0;
    for (const vector<int>& v : xvv) {
      for (const int& i : v) {
        IX[k] = six[i];
        VAL[k] = val[i];
        k++;
      }
    }

    return SparseMat<ixT, valT>(IX, P, VAL, psize);
  }
  
};

template<class ixT, class valT>
struct SparseMat {
  ixT i;
  ixT p;
  valT v;
  int pri_size;

  SparseMat(const ixT& i, const ixT& p, const valT& v, const int& pri_size)
    : i(i), p(p), v(v), pri_size(pri_size) {}

};


// SORT OF SPARSE MATRIX INDICES WITHIN COLUMNS

// adapted from http://stackoverflow.com/a/17074810/453735
template <class It>
std::vector<size_t> sort_permutation(It beg, const size_t size) {
  std::vector<size_t> p(size);
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
            [&beg](size_t i, size_t j){ return *(beg + i) < *(beg + j); });
  return p;
}

template <class It1, class It2>
void apply_permutation_in_place(It1 vec1, It2 vec2, const std::vector<size_t>& p) {
  std::vector<bool> done(p.size());
  for (size_t i = 0; i < p.size(); ++i) {
    if (done[i]) continue;
    done[i] = true;
    size_t prev_j = i;
    size_t j = p[i];
    while (i != j) {
      std::swap(vec1[prev_j], vec1[j]);
      std::swap(vec2[prev_j], vec2[j]);
      done[j] = true;
      prev_j = j;
      j = p[j];
    }
  }
}

inline void reorder(IntegerVector& IX, IntegerVector& P, NumericVector& X) {
  size_t psize = P.size();
  for (int i = 1; i < psize; i++) {
    int start = P[i - 1];
    int size = P[i] - start;
    std::vector<size_t> perm = sort_permutation(&IX[start], size);
    apply_permutation_in_place(&IX[start], &X[start], perm);
  }
}

}

#endif
