// -*- mode: c++ -*-

#ifndef DISTANCE_COMMON_HPP__
#define DISTANCE_COMMON_HPP__

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <Rcpp.h>
#include <RcppParallel.h>

using std::vector;

using namespace Rcpp;
using namespace RcppParallel;
static const double POS_INF = std::numeric_limits<double>::max();

template<class IxIt, class ValIt>
vector<double> pointwise_average(const RMatrix<double>& vecs,  IxIt ibeg, IxIt iend, ValIt ival)
{
  size_t nrows = vecs.nrow();
  vector<double> out(nrows);
  int nvecs = 0;

  while (ibeg != iend) {
    auto vec = vecs.column(*ibeg);
    auto val = (*ival);
    for (size_t vi = 0;  vi < nrows; vi++) {
      out[vi] +=  val * vec[vi];
    }
    ibeg++;
    ival++;
    nvecs++;
  }

  for (size_t i = 0; i < nrows; i++) {
    out[i] /= nvecs;
  }
  return out;
}

#endif
