#ifndef Y3_CLUSTER_DEL_SIG_TOM_HH
#define Y3_CLUSTER_DEL_SIG_TOM_HH

#include "ez.hh"
//#include "primitives.hh"
#include "interp_2d.hh"

#include <memory>
#include <cmath>

namespace y3_cluster
{
  class SIG_SUM {
  private:
    std::shared_ptr<Interp2D const> _sigma1;
    std::shared_ptr<Interp2D const> _sigma2;
    std::shared_ptr<Interp2D const> _bias;

  public:
    SIG_SUM(std::shared_ptr<Interp2D const> sigma1,
                std::shared_ptr<Interp2D const> sigma2,
                std::shared_ptr<Interp2D const> bias)
                : _sigma1(sigma1), _sigma2(sigma2), _bias(bias) {}

    using doubles = std::vector<double>;


    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    { 
      double _sig_1 = _sigma1->eval(r,lnM);
      double _sig_2 = _bias->eval(zt,lnM) * _sigma2->eval(r,zt);
      return (1.+zt)*(1.+zt)*(1.+zt)*(_sig_1+_sig_2);
    }
  };
}

#endif
