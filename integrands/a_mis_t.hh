#ifndef Y3_CLUSTER_A_MIS_T_HH
#define Y3_CLUSTER_A_MIS_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

namespace y3_cluster {

  class A_MIS_t {
  public:
    A_MIS_t() {}
    explicit A_MIS_t(cosmosis::DataBlock&) {}

    double
    operator()(double, double, double, double, double) const
    {
      // Will update with a model eventually!
      return 1.0;
    }
  };

}

#endif
