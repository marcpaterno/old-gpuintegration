#ifndef Y3_CLUSTER_A_CEN_T_HH
#define Y3_CLUSTER_A_CEN_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

namespace y3_cluster {

  class A_CEN_t {
  public:
    A_CEN_t() {}
    explicit A_CEN_t(cosmosis::DataBlock&) {}

    double
    operator()(double, double, double, double) const
    {
      // Will update with a model eventually!
      return 1.0;
    }
  };

}

#endif
