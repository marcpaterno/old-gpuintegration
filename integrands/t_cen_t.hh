#ifndef Y3_CLUSTER_T_CEN_T_HH
#define Y3_CLUSTER_T_CEN_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

namespace y3_cluster {

  class T_CEN_t {
  public:
      T_CEN_t() {}
      explicit T_CEN_t(cosmosis::DataBlock&) {}

      double
      operator()(double, double) const
      {
        // Will update with a model eventually!
        return 1.0;
      }
  };
}

#endif
