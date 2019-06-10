#ifndef Y3_CLUSTER_INTERP_1D_HH
#define Y3_CLUSTER_INTERP_1D_HH

#include "gsl/gsl_interp.h"
#include <array>
#include <cstddef>
#include <vector>

// Interp1D is used for linear interpolation in 1 dimension.
// It uses the GSL library to do the actual interpolation.
// Interp1D object allow extrapolation as well as supporting
// interpolation; no warnings or errors are given when
// extrapolating.
//
namespace y3_cluster {
  class Interp1D {
  public:
    // Interpolator created from two arrays; compiler assures they are of the
    // same length.
    template <std::size_t N>
    Interp1D(std::array<double, N> const& xs, std::array<double, N> const& ys);

    // Interpolator created from two vectors: throws std::logic_error if the
    // vectors are not of the same length.
    Interp1D(std::vector<double> const& xs, std::vector<double> const& ys);

    // Destructor must clean up allocated GSL resources.
    ~Interp1D() noexcept;

    // Interp1D objects can not be copied because the GSL resources can not
    // be copied.
    Interp1D(Interp1D const&) = delete;

    double operator()(double x) const;
    double
    eval(double x) const
    {
      return this->operator()(x);
    };

  private:
    std::vector<double> xs_;
    std::vector<double> ys_;
    gsl_interp* interp_;
  };
}

template <std::size_t N>
inline y3_cluster::Interp1D::Interp1D(std::array<double, N> const& xs,
                                      std::array<double, N> const& ys)
  : Interp1D({begin(xs), end(xs)}, {begin(ys), end(ys)})
{}

#endif
