#ifndef Y3_CLUSTER_CPP_ROFFSET_T_HH
#define Y3_CLUSTER_CPP_ROFFSET_T_HH

#include <cmath>
#include "utils/datablock_reader.hh"

namespace y3_cluster {

  class ROFFSET_t {
  public:
    explicit ROFFSET_t(double tau) : _tau(tau) {}

    explicit ROFFSET_t(cosmosis::DataBlock& sample)
        : _tau(get_datablock<double>(sample, "cluster_abundance", "roffset_tau"))
    {}

    double
    operator()(double x) const
    {
      // eq. 36
      return x / _tau / _tau * std::exp(-x / _tau);
    }

  private:
    double _tau;
  };
}

#endif
