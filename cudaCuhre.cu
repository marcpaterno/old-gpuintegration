#include <mpi.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <iomanip>  
#include "quad/quad.h"
#include "function.cu"

#include "quad/util/cudaUtil.h"

#include "quad/GPUquad/GPUquad.cu"
#include "lc_lt_t.h"

/*
//cubacpp requires c++17 to compile
#include "cubacpp/cubacpp.hh"
#include "cubacpp/array.hh"
#include "cubacpp/cuhre.hh"
#include "cubacpp/integration_volume.hh"
#include "cubacpp/gsl.hh"
#include "lc_lt_t.hh"
*/

//#include "CuhreInteface.h"
#include <chrono>

using namespace quad;

int main(int argc, char **argv){
  
  TYPE epsrel = 1e-12;
  TYPE epsabs = 1e-12;
  int verbose = 0;
  
  //Num Devices 
  int numDevices = 1;
  //auto start = std::chrono::steady_clock::now();	
   
  TYPE integral = 0, error = 0;
  size_t nregions = 0, neval = 0;
  int dim = 3;

  GPUcuhre<TYPE> *cuhre = new GPUcuhre<TYPE>(argc, argv, dim, 0, verbose, numDevices);
  
  //GPUcuhre instantiation causes the device to reset. Instantiating LC_LT_t allocates device memory
  //so it must occur after the GPUcuhre instantiation.
  LC_LT_t<IntegralGPU> t;
  
  float highs[3] = {2, .3, 1};
  float lows[3] =  {1, .1, 0};
  Volume vol;
  vol.Initialize(lows, highs, 3);
  
  int errorFlag = cuhre->integrate(epsrel, epsabs, integral, error, nregions, neval, &t);
    
  printf("%d\t%e\t%.10lf\t%.10f\t%ld\t%ld\t%d\n", dim, epsrel, integral, error, nregions, neval, errorFlag);

  delete cuhre;
    
  return 0;
}
