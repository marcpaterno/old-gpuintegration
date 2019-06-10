#ifndef Y3_CLUSTER_SIGMA_HH
#define Y3_CLUSTER_SIGMA_HH

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "utils/bessel_polynomial_integral.hh"
#include "utils/datablock_reader.hh"
#include "utils/integration_range.hh"
#include "utils/interp_1d.hh"
#include "utils/interp_2d.hh"
#include "utils/primitives.hh"
#include "utils/read_vector.hh"

#include "omega_z_sdss.hh"

namespace y3_cluster {
  // Coefficients of survey mask, in spherical harmonics
  template<typename OMEGA_Z>
  double
  survey_mask_kay(const OMEGA_Z& omega_z, unsigned l, double z)
  {
    if (l == 0)
        return 0.5 / std::sqrt(pi());
    else {
      const double oz = omega_z(z),
                   cos_theta = 1 - oz / 2 / pi(),
                   frac = (std::legendre(l - 1, cos_theta) - std::legendre(l + 1, cos_theta))
                          / oz;
      return std::sqrt(pi() / (2*l + 1)) * frac;
    }
  }

  namespace {
    // Create a grid of K_l(z_i), so: compute_ks()[i][j] = K_j(z_i)
    template<typename OMEGA_Z>
    static std::vector<std::vector<double>>
    compute_ks(const OMEGA_Z& omega_z, std::size_t maxl, const std::vector<IntegrationRange>& zirs)
    {
      std::vector<std::vector<double>> out;
      for (const auto zbin : zirs) {
        std::vector<double> this_bin;
        // Use the midpoint of the bin
        for (auto l = 0u; l < maxl; l++)
          this_bin.push_back(survey_mask_kay(omega_z, l, zbin.transform(0.5)));
        out.push_back(this_bin);
      }
      return out;
    }
  }

  /* The Sample Variance covariance is computed by:
   *
   *  C_{ij} = Nb_i * Nb_j * sigma_{ij}^2
   *
   * The class in this file computes a matrix of the last term - sigma_{ij}^2
   */
  class SampleVariance_t {
    std::vector<IntegrationRange> z_ranges;
    std::shared_ptr<Interp2D const> matter_power_lin;
    std::shared_ptr<Interp1D const> dcom;

    // Number of Bessel functions to sum over, i.e. l = 0, 1, ..., maxl
    std::size_t maxl;
    std::vector<std::vector<double>> ks;

    // Reduced hubble parameter `h`
    double hubble;

  public:
    SampleVariance_t() = delete;

    template<typename OMEGA_Z>
    SampleVariance_t(const OMEGA_Z& omega_z, const std::vector<IntegrationRange>& ir, double h)
      : z_ranges(ir)
      , matter_power_lin(std::make_shared<Interp2D const>(
            read_vector("matter_power_lin/k_h.txt"),
            read_vector("matter_power_lin/z.txt"),
            read_vector("matter_power_lin/p_k.txt")))
      , dcom(std::make_shared<Interp1D const>(
            read_vector("distances/z.txt"),
            read_vector("distances/d_m.txt")))
      , maxl(90)
      , ks(compute_ks(omega_z, maxl, ir))
      , hubble(h)
      {}

    // TODO point to arbitrary matter_power_lin section name
    template<typename OMEGA_Z>
    explicit SampleVariance_t(cosmosis::DataBlock& sample, const OMEGA_Z& omega_z, const std::vector<IntegrationRange>& z_ranges)
      : z_ranges(z_ranges)
      , matter_power_lin(std::make_shared<Interp2D const>(
            get_datablock<std::vector<double>>(sample, "matter_power_lin_cdm_baryon", "k_h"),
            get_datablock<std::vector<double>>(sample, "matter_power_lin_cdm_baryon", "z"),
            get_datablock<cosmosis::ndarray<double>>(sample, "matter_power_lin_cdm_baryon", "p_k")))
      , dcom(std::make_shared<Interp1D const>(
            get_datablock<std::vector<double>>(sample, "distances", "z"),
            get_datablock<std::vector<double>>(sample, "distances", "d_m")))
      , maxl(get_datablock<int>(sample, "cluster_abundance", "smp_var_maxl"))
      , ks(compute_ks(omega_z, maxl, z_ranges))
      , hubble(get_datablock<double>(sample, "cosmological_parameters", "h0"))
      {}

    /* Computes:
     *
     *    sum_{l=0}^{maxl} K_l(i) * K_l(j) * R_i(k) * R_j(k)
     */
    double
    compute_sum_over_bessels(double k, unsigned i, unsigned j) const;

    double
    manual_compute_sum_over_bessels(double k, unsigned i, unsigned j) const;

    // Compute the full sigma_{ij} matrix
    std::vector<std::vector<double>> compute() const;
  };
}

#endif
