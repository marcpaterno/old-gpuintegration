#include <cubacpp/cubacpp.hh>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include "sample_variance.hh"

using y3_cluster::IntegrationRange;

/* Computes:
 *
 *    sum_{l=0}^{maxl} K_l(i) * K_l(j) * R_i(k) * R_j(k)
 */
double
y3_cluster::SampleVariance_t::compute_sum_over_bessels(double k,
                                                       unsigned i,
                                                       unsigned j) const
{
  const double ri_min = dcom->eval(z_ranges[i].transform(0)) * hubble,
               ri_max = dcom->eval(z_ranges[i].transform(1)) * hubble,
               rj_min = dcom->eval(z_ranges[j].transform(0)) * hubble,
               rj_max = dcom->eval(z_ranges[j].transform(1)) * hubble,
               i_normalization =
                 (ri_max * ri_max * ri_max - ri_min * ri_min * ri_min) / 3.0,
               j_normalization =
                 (rj_max * rj_max * rj_max - rj_min * rj_min * rj_min) / 3.0;

  const std::vector<double> bessels_i = bessel_polynomial_integrals(
                              2, maxl, k, ri_min, ri_max),
                            bessels_j = (i == j) ?
                                          bessels_i :
                                          bessel_polynomial_integrals(
                                            2, maxl, k, ri_min, ri_max);

  double sum = 0;
  for (auto l = 0u; l < maxl; l++) {
    sum += bessels_i[l] * bessels_j[l] / i_normalization / j_normalization *
           ks[i][l] * ks[j][l];
  }
  return sum;
}

double
y3_cluster::SampleVariance_t::manual_compute_sum_over_bessels(double k,
                                                              unsigned i,
                                                              unsigned j) const
{
  const double ri_min = dcom->eval(z_ranges[i].transform(0)) * hubble,
               ri_max = dcom->eval(z_ranges[i].transform(1)) * hubble,
               rj_min = dcom->eval(z_ranges[j].transform(0)) * hubble,
               rj_max = dcom->eval(z_ranges[j].transform(1)) * hubble,
               i_normalization =
                 (ri_max * ri_max * ri_max - ri_min * ri_min * ri_min) / 3.0,
               j_normalization =
                 (rj_max * rj_max * rj_max - rj_min * rj_min * rj_min) / 3.0;

  cubacpp::Cuhre integrator;
  integrator.maxeval = 999999999;
  auto integrate_bessel = [&](
                            unsigned l, double k, double r_min, double r_max) {
    y3_cluster::IntegrationRange r_ir{r_min, r_max};
    const auto result = integrator.integrate(
      [&](double r_scaled, double ignored [[maybe_unused]]) {
        const double r = r_ir.transform(r_scaled);
        return r * r * gsl_sf_bessel_jl(l, k * r) * r_ir.jacobian();
      },
      1e-2,
      1e-18);

    if (result.status != 0)
      throw std::runtime_error("Big boom! (bessel integral did not converge)");

    return result.value;
  };

  double sum = 0;
  for (auto l = 0u; l < maxl; l++) {
    sum += integrate_bessel(l, k, ri_min, ri_max) *
           integrate_bessel(l, k, rj_min, rj_max) / i_normalization /
           j_normalization * ks[i][l] * ks[j][l];
  }
  return sum;
}

std::vector<std::vector<double>>
y3_cluster::SampleVariance_t::compute() const
{
  std::vector<std::vector<double>> sigma_mat(z_ranges.size());
  for (auto& v : sigma_mat)
    v.resize(z_ranges.size());

  // Compute matrix of sigma^{SampVar}_{ij}
  cubacpp::Cuhre integrator;
  integrator.maxeval = 999999999;

  // Using Matteo's range for now
  IntegrationRange lnk_ir{std::log(0.0001), std::log(0.80)};
  for (auto i = 0u; i < z_ranges.size(); i++) {
    for (auto j = 0u; j < z_ranges.size(); j++) {
      if (i <= j) {
        /* Integrate \sigma_{ij} =
         *      ((dV/d\Omega)|_{\Delta zi} (dV/d\Omega)|_{\Delta zi})^{-1}
         *       2/pi \int dk^3 \int_{\Delta z_i} dz_i \int_{\Delta z_j} dz_j
         *       \sqrt{P_{lin}(k, z_i) P_{lin}(k, z_j)} D(z_i)^2 D(z_j)^2
         *       \sum_l (-1}^l j_l(k D(z_i)) j_l(k D(z_j)) K_l(z_i) K_l(z_j)
         */

        // Have ignored parameter so we can use Cuhre - which requires >=2
        // integration dimensions But, Cuhre is much, much faster at this than
        // Vegas, Suave
        auto result = integrator.integrate(
          [&](const double lnk_scaled, const double ignored [[maybe_unused]]) {
            const auto lnk = lnk_ir.transform(lnk_scaled),
                       zi_mid = z_ranges[i].transform(0.5),
                       zj_mid = z_ranges[j].transform(0.5), k = std::exp(lnk);

            // TODO Matteo computes the sum-over-bessels ahead of
            // time on a grid in ln(k) and interpolates over it -
            // we should do the same
            const auto matter_power =
                         std::sqrt(matter_power_lin->eval(k, zi_mid) *
                                   matter_power_lin->eval(k, zj_mid)),
                       sum = compute_sum_over_bessels(k, i, j);

            // k^3 due to integration over log
            return 2 / pi() * k * k * k * sum * matter_power *
                   lnk_ir.jacobian();
          },
          1e-2,
          1e-18);

        if (result.status != 0)
          throw std::runtime_error(
            "Integration of dvdodz failed in SampleVariance");

        sigma_mat[i][j] = result.value;
      } else {
        // Already computed (sigma_{ij} = sigma_{ji})
        sigma_mat[i][j] = sigma_mat[j][i];
      }
    }
  }
  return sigma_mat;
}
