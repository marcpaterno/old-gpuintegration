//#include <thrust/device_vector.h>
//#include <thrust/host_vector.h> 
#include "lc_lt_t.h"
#include "quad/GPUquad/GPUquad.cu"
#include <chrono>

using namespace quad;

class test {
public:
  __device__ double
  operator()(double xx[], int dim)
  {
    double f = 0;
    double sum = 0;
    int N = 1;
    for (N = 1; N <= dim; ++N)
      sum = sum - cos(10.0 * xx[N - 1]) / 0.054402111088937;

    f = sum / 2.0;
    return f;
  }
};

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

  // GPUcuhre<TYPE> *cuhre = new GPUcuhre<TYPE>(argc, argv, dim, 0, verbose,
  // numDevices);
  GPUcuhre<TYPE> cuhre(argc, argv, dim, 0, verbose, numDevices);
  // GPUcuhre instantiation causes the device to reset. Instantiating LC_LT_t
  // allocates device memory so it must occur after the GPUcuhre instantiation.
  LC_LT_t<IntegralGPU> t;

  float highs[3] = {2, .3, 1};
  float lows[3] = {1, .1, 0};
  Volume vol(lows, highs, 3);

  // int errorFlag = cuhre->integrate(epsrel, epsabs, integral, error, nregions,
  // neval, &t, &vol); int errorFlag = cuhre->integrate(epsrel, epsabs, integral,
  // error, nregions, neval, &t);
  int errorFlag =
    cuhre.integrate(epsrel, epsabs, integral, error, nregions, neval, &t, &vol);
  printf("%d\t%e\t%.10lf\t%.10f\t%ld\t%ld\t%d\n",
         dim,
         epsrel,
         integral,
         error,
         nregions,
         neval,
         errorFlag);

  printf("=================\n");
  integral = 0;
  error = 0;
  nregions = 0;
  neval = 0;
  errorFlag =
    cuhre.integrate(1e-7, epsabs, integral, error, nregions, neval, &t, &vol);
  // int errorFlag = cuhre->integrate(epsrel, epsabs, integral, error, nregions,
  // neval, &t);

  printf("%d\t%e\t%.10lf\t%.10f\t%ld\t%ld\t%d\n",
         dim,
         1e-7,
         integral,
         error,
         nregions,
         neval,
         errorFlag);

  // delete cuhre;
  return 0;
}
