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

#endif
