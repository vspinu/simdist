#ifndef DISTANCE_NORMS_HPP__
#define DISTANCE_NORMS_HPP__

namespace simdist {

enum NormType {NONORM, L1, L2};

template <NormType NT>
struct Norm {
  Norm() {}
  virtual ~Norm() {}

  template <class ValIt>
  inline void operator()(ValIt xbeg, ValIt xend) const {
    throw std::invalid_argument("Unsupported NormType");
  };  
};

template <>
struct Norm<L2> {
  template <class FwIt> 
  inline void operator()(FwIt xbeg, FwIt xend) const {
    double sum = 0.0;
    FwIt xbeg2 = xbeg;
    while (xbeg != xend) {
      sum += (*xbeg) * (*xbeg);
      xbeg++;
    }
    sum = std::sqrt(sum);
    if (sum != 0.0) {
      while (xbeg2 != xend) {
        *xbeg2 /= sum;
        xbeg2++;
      }
    }
  }
};

}

#endif
