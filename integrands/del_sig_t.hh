#ifndef Y3_CLUSTER_DEL_SIG_T_HH
#define Y3_CLUSTER_DEL_SIG_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "ez.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include "utils/primitives.hh"

#include <cmath>
#include <memory>

namespace y3_cluster {
  class DEL_SIG_t {
  private:
    std::shared_ptr<Interp2D const> _dsigma1;
    std::shared_ptr<Interp2D const> _dsigma2;
    std::shared_ptr<Interp2D const> _bias;

  public:
    DEL_SIG_t(std::shared_ptr<Interp2D const> dsigma1,
              std::shared_ptr<Interp2D const> dsigma2,
              std::shared_ptr<Interp2D const> bias)
      : _dsigma1(dsigma1), _dsigma2(dsigma2), _bias(bias)
    {}

    using doubles = std::vector<double>;

    // TODO: This needs to be reading cosmosis datablock parameters
    explicit DEL_SIG_t(cosmosis::DataBlock& sample)
      : _dsigma1(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_x1"),
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_y1"),
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_z1")))
      , _dsigma2(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_x2"),
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_y2"),
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_z2")))
      , _bias(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_x3"),
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_y3"),
          get_datablock<doubles>(sample,
                                 "cluster_abundance",
                                 "del_sig_params_z3")))
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double del_sig_1 = _dsigma1->eval(r, lnM);
      double del_sig_2 = _bias->eval(zt, lnM) * _dsigma2->eval(r, zt);
      if (del_sig_1 >= del_sig_2) {
        return (1. + zt) * (1. + zt) * (1. + zt) * del_sig_1;
      } else {
        return (1. + zt) * (1. + zt) * (1. + zt) * del_sig_2;
      }
    }
  };
}

#endif
