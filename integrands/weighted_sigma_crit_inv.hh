#ifndef Y3_CLUSTER_WEIGHTED_SIG_CRIT_INV
#define Y3_CLUSTER_WEIGHTED_SIG_CRIT_INV

#include <exception>
#include <memory>
#include <vector>

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "cubacpp/gsl.hh"
#include "sigma_crit_inverse_t.hh"
#include "utils/datablock_reader.hh"
#include "utils/fits_loader.hh"
#include "utils/interp_1d.hh"
#include "utils/primitives.hh"

/* Computes, for each lensing source distribution P_i(\zs), the integral
 *
 * \Sigma_{w, crit}^{-1}(\zt)
 *      = \int_{\zt}^\inf d(\zs) \Sigma_{crit}^{-1}(\zt, \zs) P(\zs)
 *      = 4 \pi G D(\zt) / (c * c) \int_{\zt}^{\inf} d(\zs) (1 - \frac{(1 + \zt)
 * D(\zt)}{(1 + \zs) D(\zs)} P(\zs) = 4 \pi G D(\zt) / (c * c)
 * [\int_{\zt}^{\inf} P(\zs) d(\zs)
 *                                 - (1 + \zt) D(\zt) \int_{\zt}^{\inf} P(\zs)/
 * ((1 + \zs) * D(\zs)) d(\zs) ]
 *
 * (Plug that in LaTeX - hard to read.)
 */
namespace y3_cluster {
  namespace {
    std::vector<std::shared_ptr<Interp1D const>>
    compute_weighted_sigma_crits(std::vector<double> distance_zs,
                                 std::vector<double> da,
                                 std::vector<double> src_zs,
                                 std::vector<std::vector<double>> pzsrcs)
    {
      // Distance metric interpolator
      sigma_crit_inv sci(std::make_shared<Interp1D const>(distance_zs, da));
      // Integrator - adaptive quadrature should be good
      cubacpp::QAG igrtr(0.0, 1.0, GSL_INTEG_GAUSS61, 1000);

      std::vector<std::shared_ptr<Interp1D const>> out;
      for (auto srci = 0u; srci < pzsrcs.size(); srci++) {
        const Interp1D pzsrc(src_zs, pzsrcs[srci]);

        std::vector<double> wgted_sig_crits;
        for (auto zt : src_zs) {
          if (pzsrc.eval(zt) == 0.0) {
            wgted_sig_crits.push_back(0.0);
            continue;
          }

          auto this_igrnd = igrtr.with_range(zt, src_zs[src_zs.size() - 1]);
          auto const wgted_sci = this_igrnd.integrate(
            [&](double zs) { return sci(zt, zs) * pzsrc.eval(zs); },
            1e-6,
            1e-12);

          if (wgted_sci.status != 0)
            throw std::runtime_error(
              "weighted sigma_crit_inv did not converge");

          wgted_sig_crits.push_back(wgted_sci.value);
        }

        out.push_back(
          std::make_shared<Interp1D const>(src_zs, wgted_sig_crits));
      }

      return out;
    }
  }

  class weighted_sigma_crit_inv {
    std::vector<double> zoffsets;
    std::vector<std::shared_ptr<Interp1D const>> pzsrcs;

  public:
    weighted_sigma_crit_inv(std::vector<double> offsets,
                            std::vector<double> dist_zs,
                            std::vector<double> da,
                            std::vector<double> src_zs,
                            std::vector<std::vector<double>> pzs)
      : zoffsets(offsets)
      , pzsrcs(compute_weighted_sigma_crits(dist_zs, da, src_zs, pzs))
    {}

    // The pzsource distrs are a separate argument so that the dists
    // only need to be read from a fits file once - unnecessary file
    // io will really hurt when running on big clusters!
    explicit weighted_sigma_crit_inv(cosmosis::DataBlock& sample,
                                     std::vector<double> src_zs,
                                     std::vector<std::vector<double>> pzsources)
      // TODO - the number of source bins should be flexible
      : weighted_sigma_crit_inv(
          {get_datablock<double>(sample,
                                 "cluster_abundance",
                                 "pzsource_offset_1"),
           get_datablock<double>(sample,
                                 "cluster_abundance",
                                 "pzsource_offset_2"),
           get_datablock<double>(sample,
                                 "cluster_abundance",
                                 "pzsource_offset_3"),
           get_datablock<double>(sample,
                                 "cluster_abundance",
                                 "pzsource_offset_4")},
          get_datablock<std::vector<double>>(sample, "distances", "z"),
          get_datablock<std::vector<double>>(sample, "distances", "d_a"),
          src_zs,
          pzsources)
    {}

    std::size_t
    nsources() const
    {
      return pzsrcs.size();
    }

    double
    operator()(std::size_t bin_no, double zt) const
    {
      return pzsrcs[bin_no]->eval(zt + zoffsets[bin_no]);
    }
  };
}

#endif
