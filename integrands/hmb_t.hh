#ifndef Y3_CLUSTER_CPP_BMZ_T_HH
#define Y3_CLUSTER_CPP_BMZ_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/ndarray.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include "utils/read_vector.hh"

#include <memory>

namespace y3_cluster {
  // Bias for the Halo Mass Function (Halo Mass Bias)
  class HMB_t {
    std::shared_ptr<Interp2D const> _tinker_bias;

  public:
    HMB_t()
      : _tinker_bias(std::make_shared<Interp2D const>(
          read_vector("tinker_bias_function/ln_mass_h.txt"),
          read_vector("tinker_bias_function/z.txt"),
          read_vector("tinker_bias_function/bias.txt")))
    {}

    explicit HMB_t(cosmosis::DataBlock& sample)
      : _tinker_bias(std::make_shared<Interp2D const>(
          get_datablock<std::vector<double>>(sample,
                                             "tinker_bias_function",
                                             "ln_mass_h"),
          get_datablock<std::vector<double>>(sample,
                                             "tinker_bias_function",
                                             "z"),
          get_datablock<cosmosis::ndarray<double>>(sample,
                                                   "tinker_bias_function",
                                                   "bias")))
    {}

    double
    operator()(double lnM, double zt) const
    {
      return _tinker_bias->eval(lnM, zt);
    }
  };
}

#endif
