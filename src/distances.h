#ifndef DISTANCE_DISTANCES_HPP__
#define DISTANCE_DISTANCES_HPP__

#include <unordered_map>

namespace simdist {

enum DistType {NODIST, COSINE, EUCLIDEAN};
enum AggrType {NOAGGR, CENTROID, SEM_MIN_MAX, SEM_MIN_SUM};

template <DistType DT>
struct Dist {
  Dist() {}
  virtual ~Dist() {}
  template <class ValIt>
  inline double operator()(ValIt xbeg, ValIt xend, ValIt ybeg) const {
    throw std::invalid_argument("Unsupported DistType");
  };  
};

template <>
struct Dist<COSINE> {

  template <class ValIt>
  inline double operator()(ValIt xbeg, ValIt xend, ValIt ybeg) const {
    return 1.0 - std::inner_product(xbeg, xend, ybeg, 0.0);
  }

  template <class IxIt, class ValIt>
  inline double operator()(IxIt xix_beg, IxIt xix_end, ValIt xval,
                           IxIt yix_beg, IxIt yix_end, ValIt yval) const {

    if (xix_beg == xix_end && yix_beg == yix_end)
      return 1.0;
    
    double out = 0.0;
    while (xix_beg != xix_end && yix_beg != yix_end) {
      /* printf("ENTRY: XVAL:%f  XVAL:%f xix:%d yix:%d\n", *xval, *yval, *xix_beg, *yix_beg); */
      if (*xix_beg == *yix_beg) {
        /* printf("=: xval:%f  yval:%f\n", *xval, *yval); */
        out += (*xval) * (*yval);
        xix_beg++; xval++;
        yix_beg++; yval++;
      }
      while (yix_beg != yix_end && *xix_beg > *yix_beg) {
        /* printf(">: xval:%f  yval:%f\n", *xval, *yval); */
        yix_beg++; yval++;
      }
      while (xix_beg != xix_end && *xix_beg < *yix_beg) {
        /* printf("<: xval:%f  yval:%f\n", *xval, *yval); */
        xix_beg++; xval++;
      }
    }
    return 1.0 - out;
  }
};

template <>
struct Dist<EUCLIDEAN> {

  template <class ValIt>
  inline double operator()(ValIt xbeg, ValIt xend, ValIt ybeg) const {
    double out = 0.0;
    while (xbeg != xend) {
      double diff = *xbeg - *ybeg;
      out += diff * diff;
      xbeg++;
      ybeg++;
    }
    return std::sqrt(out);
  }
  
  template <class IxIt, class ValIt>
  inline double operator()(IxIt xix_beg, IxIt xix_end, ValIt xval,
                           IxIt yix_beg, IxIt yix_end, ValIt yval) const {
    double out = 0.0;
    while (xix_beg != xix_end && yix_beg != yix_end) {
      /* printf("ENTRY: XVAL:%f  XVAL:%f xix:%d yix:%d\n", *xval, *yval, *xix_beg, *yix_beg); */
      if (*xix_beg == *yix_beg) {
        /* printf("=: xval:%f  yval:%f\n", *xval, *yval); */
        double diff = *xval - *yval;
        out += diff * diff;
        xix_beg++; xval++;
        yix_beg++; yval++;
      }
      while (yix_beg != yix_end && *xix_beg > *yix_beg) {
        /* printf(">: xval:%f  yval:%f\n", *xval, *yval); */
        out += (*yval) * (*yval);
        yix_beg++; yval++;
      }
      while (xix_beg != xix_end && *xix_beg < *yix_beg) {
        /* printf("<: xval:%f  yval:%f\n", *xval, *yval); */
        out += (*xval) * (*xval);
        xix_beg++; xval++;
      }
    }

    while (yix_beg != yix_end) {
      out += (*yval) * (*yval);
      yix_beg++; yval++;
    }
    while (xix_beg != xix_end) {
      out += (*xval) * (*xval);
      xix_beg++; xval++;
    }
    
    return std::sqrt(out);
  }

};

} // end simdist namespace

#endif
