#ifndef Y3_CLUSTER_T_MIS_T_HH
#define Y3_CLUSTER_T_MIS_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

namespace y3_cluster {

  class T_MIS_t {
  public:
      T_MIS_t() {}
      explicit T_MIS_t(cosmosis::DataBlock&) {}

      double
      operator()(double) const
      {
        // Will update with a model eventually!
        return 1.0;
      }
  };
}

#endif
