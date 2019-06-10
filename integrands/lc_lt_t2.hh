#ifndef Y3_CLUSTER_CPP_LC_LT_T2_HH
#define Y3_CLUSTER_CPP_LC_LT_T2_HH

#include "utils/primitives.hh"

namespace y3_cluster {

  struct LC_LT_t2 {
    double
    operator()(double lc, double lt, double /* zt */) const
    {
      return y3_cluster::gaussian(lc - lt, 0.0, 10.0);
    }
  };
}

#endif
