#ifndef Y3_CLUSTER_LC_LT_T_HH
#define Y3_CLUSTER_LC_LT_T_HH

//#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "interp_2d.hh"
#include "primitives.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace y3_cluster {
  struct LC_LT_t {

    static Interp2D const tau_interp;
    static Interp2D const mu_interp;
    static Interp2D const sigma_interp;
    static Interp2D const fmsk_interp;
    static Interp2D const fprj_interp;

    explicit LC_LT_t() {}
    //LC_LT_t() {}

    double
    operator()(double lc, double lt, double zt) const
    {
      //lc = 0.281076;
      //lt = 1.977469;
      //zt = 0.156215;
      const auto tau = tau_interp(lt, zt);
      const auto mu = mu_interp(lt, zt);
      const auto sigma = sigma_interp(lt, zt);
      const auto fmsk = fmsk_interp(lt, zt);
      const auto fprj = fprj_interp(lt, zt);
      
      //const auto mu = 1;
      //const auto sigma =.5;
      //const auto fmsk = .6;
      //const auto fprj =.6;
	//printf("lt:%f zt:%f tau:%f\n", 2, .12, tau_interp(2, .12));
	//printf("lt:%f zt:%f mu:%f\n", 2, .12, mu_interp(2, .12));
	
	
	//printf("lt:%f zt:%f fprj:%f\n", 2, .12, fprj_interp(2, .12));
	//auto t = tau_interp(2, .12);
      //printf("Result for 2 .12: %f\n", t);
      //t = tau_interp(1, .17);
      //printf("Result for 1 .17: %f\n", t);
      
      //t = tau_interp(1, 1);
      //printf("Result for 1 1: %f\n", t);

      const auto exptau =
        std::exp(tau * (2.0 * mu + tau * sigma * sigma - 2.0 * lc) / 2.0);
      const auto root_two_sigma = std::sqrt(2.0) * sigma;
      const auto mu_tau_sig_sqr = mu + tau * sigma * sigma;
      //printf("root_two_sigma:%f\n", root_two_sigma);
      //printf("exptau:%f\n", exptau);
      //printf("mu_tau_sig_sqr:%f\n", mu_tau_sig_sqr);
     
      // Helper function for common pattern
      const auto erfc_scaled = [root_two_sigma](double a, double b) {
          return std::erfc((a - b) / root_two_sigma);
      };
      //printf("erfc1:%f\n", erfc_scaled(lc, mu));
      //printf("gauss1:%f\n",  y3_cluster::gaussian(lc, mu, sigma));
      // eq. (33)
      auto f = (1.0 - fmsk) * (1.0 - fprj) * y3_cluster::gaussian(lc, mu, sigma) +
             0.5 * ((1.0 - fmsk) * fprj * tau + fmsk * fprj / lt) * exptau * erfc_scaled(mu_tau_sig_sqr, lc) +
             0.5 * fmsk / lt *(erfc_scaled(lc, mu) - erfc_scaled(lc + lt, mu)) -
             0.5 * fmsk * fprj / lt *(std::exp(-tau * lt) * exptau *erfc_scaled(mu_tau_sig_sqr, lc + lt));
      //printf("feval:%f\n", f);
      //printf("lt: %f zt:%f lc:%f tau:%f feval:%f\n",lt, zt, lc,  tau, f);
      return f;
    }
  };
}

#endif
