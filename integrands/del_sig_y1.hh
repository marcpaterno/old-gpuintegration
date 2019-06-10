#ifndef Y3_CLUSTER_DEL_SIG_CEN_Y1_HH
#define Y3_CLUSTER_DEL_SIG_CEN_Y1_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include <cmath>

namespace y3_cluster

{
  /* NFW Profile , in h*M_solar*Mpc^-2 */
  struct DEL_SIG_y1 {

    DEL_SIG_y1() {}
    DEL_SIG_y1(const cosmosis::DataBlock&) {}

    double
    operator()(double, double lnM, double) const
    {
      return std::exp(lnM);
    }
  };
}

#endif
