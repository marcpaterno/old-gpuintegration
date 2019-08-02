#ifndef Y3_CLUSTER_PZSOURCE_HH
#define Y3_CLUSTER_PZSOURCE_HH

#include <stdexcept>

#include "utils/datablock_reader.hh"
#include "utils/primitives.hh"

namespace y3_cluster {
  class PZSOURCE_t {
  public:
    // TODO take in pzsource distributions
    explicit PZSOURCE_t(
      std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources,
      std::vector<double> zoffsets)
      : pzsources(pzsources), zoffsets(zoffsets)
    {}

    explicit PZSOURCE_t(
      cosmosis::DataBlock& sample,
      std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources)
      : pzsources(pzsources)
      , zoffsets{get_datablock<double>(sample,
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
                                       "pzsource_offset_4")}
    {
      if (pzsources.size() != zoffsets.size())
        throw std::runtime_error(
          "number of zsource distributions and number of offsets do not match");
    }

    double
    operator()(std::size_t bin_no, double zs) const
    {
      return pzsources[bin_no]->eval(zs + zoffsets[bin_no]);
    }

    std::size_t
    nsources() const
    {
      return pzsources.size();
    }

  private:
    std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources;
    std::vector<double> zoffsets;
  };
}

#endif
