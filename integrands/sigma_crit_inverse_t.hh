#ifndef Y3_CLUSTER_SIG_CRIT_INV
#define Y3_CLUSTER_SIG_CRIT_INV

#include <exception>
#include <memory>
#include <vector>

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_1d.hh"
#include "utils/primitives.hh"

namespace y3_cluster {

  class sigma_crit_inv {
  public:
    sigma_crit_inv(std::shared_ptr<Interp1D const> da)
      : _da(da)
      // Constants are in primitives.hh
      // Units: Mpc/s (NO h)
      , _c(c())
      // Units: Mpc^3 / M_sol / s^2
      , _G(G())
    {}

    explicit sigma_crit_inv(cosmosis::DataBlock& sample)
      : sigma_crit_inv(std::make_shared<y3_cluster::Interp1D const>(
                     get_datablock<std::vector<double>>(sample, "distances", "z"),
                     get_datablock<std::vector<double>>(sample, "distances", "d_a")))
    {}

    double
    operator()(double zt, double zs) const
    {
      double _sig_crit_inv = 0;
      if (zs > zt) {
        // Units in Mpc, M_sol, s
        double const da_zt = _da->eval(zt), // da_z needs to be in Mpc
                     da_zs = _da->eval(zs), // da_z needs to be in Mpc
                     da_zt_zs = da_zs - (1.0+zt)/(1.0+zs) *da_zt; // da_z needs to be in Mpc
        _sig_crit_inv = 4.0*pi()*G()/c()/c() * da_zt * da_zt_zs / da_zs;
      }
      // Mpc^2 / M_sol
      return _sig_crit_inv;
    }

  private:
    std::shared_ptr<Interp1D const> _da;
    double _c;
    double _G;
  };
}

#endif
