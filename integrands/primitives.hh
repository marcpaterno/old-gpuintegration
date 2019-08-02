#ifndef Y3_CLUSTER_CPP_PRIMITIVES_HH
#define Y3_CLUSTER_CPP_PRIMITIVES_HH

// include <polynomial.hh>

#include <cmath>
#include <fstream>
#include <iostream>
// primitives.hh contains a few commonly-used mathematical primitives.

namespace y3_cluster {

  inline double
  pi()
  {
    return 4. * std::atan(1.0);
  };

  inline double
  invsqrt2pi()
  {
    return 1. / std::sqrt(2. * pi());
  };

  inline double
  gaussian(double x, double mu, double sigma)
  {
    double const z = (x - mu) / sigma;
    return std::exp(-z * z / 2.) * 0.3989422804014327 / sigma;
  }

  namespace {
    // Tail recursive helper for `integer_pow`
    constexpr double
    do_integer_pow(const double accumulator, const double n, const unsigned pow)
    {
      if (pow == 0)
        return accumulator;
      if ((pow % 2) == 0)
        return do_integer_pow(accumulator, n * n, pow / 2);
      return do_integer_pow(accumulator * n, n, pow - 1);
    }
  }

  // In C++ >= 11, the std::pow does not optimize for integers :/
  constexpr double
  integer_pow(double n, int pow)
  {
    if (pow == 0)
      return 1;
    if (pow < 0)
      return do_integer_pow(1, 1.0 / n, 0 - pow);
    return do_integer_pow(1, n, pow);
  }
}

#endif
