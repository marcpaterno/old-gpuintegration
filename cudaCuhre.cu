#include <mpi.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "function.cu"
#include "quad/quad.h"
#include <iomanip>

#include "quad/util/cudaUtil.h"

#include "lc_lt_t.h"
#include "quad/GPUquad/GPUquad.cu"

/*
//cubacpp requires c++17 to compile
#include "cubacpp/array.hh"
#include "cubacpp/cubacpp.hh"
#include "cubacpp/cuhre.hh"
#include "cubacpp/gsl.hh"
#include "cubacpp/integration_volume.hh"
#include "lc_lt_t.hh"
*/

//#include "CuhreInteface.h"
#include <chrono>

using namespace quad;

int
main(int argc, char** argv)
{

  TYPE epsrel = 1e-12;
  TYPE epsabs = 1e-12;
  int verbose = 0;

  // Num Devices
  int numDevices = 1;
  // auto start = std::chrono::steady_clock::now();

  TYPE integral = 0, error = 0;
  size_t nregions = 0, neval = 0;
  int dim = 3;

  GPUcuhre<TYPE> cuhre(argc, argv, dim, 0, verbose, numDevices);

  // GPUcuhre instantiation causes the device to reset. Instantiating LC_LT_t
  // allocates device memory so it must occur after the GPUcuhre instantiation.
  LC_LT_t<IntegralGPU> t;

  int errorFlag =
    cuhre.integrate(epsrel, epsabs, integral, error, nregions, neval, &t);
  printf("%d\t%e\t%.10lf\t%.10f\t%ld\t%ld\t%d\n",
         dim,
         epsrel,
         integral,
         error,
         nregions,
         neval,
         errorFlag);

  int errorFlag2 = cuhre.integrate(
    epsrel / 10, epsabs / 10, integral, error, nregions, neval, &t);
  printf("%d\t%e\t%.10lf\t%.10f\t%ld\t%ld\t%d\n",
         dim,
         epsrel,
         integral,
         error,
         nregions,
         neval,
         errorFlag2);
  return 0;
}
