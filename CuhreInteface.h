#ifndef CUHREINTERFACE_H
#define CUHREINTERFACE_H

#include "quad/GPUquad/GPUQuadInterp2D.cu"
//#include "cubacpp/cubacpp.hh"
#include "../../gpuintegration/integrands/interp_2d.hh"

struct IntegralGPU {
  typedef quad::Interp2D Interp2D;
  typedef quad::GPUcuhre<double> Cuhre;
};

struct IntegralCPU {
  typedef y3_cluster::Interp2D Interp2D;
  // typedef cubacpp::Cuhre Cuhre; //c++17 required
};

__device__ double
gaussian(double x, double mu, double sigma)
{
  double z = (x - mu) / sigma;
  return exp(-z * z / 2.) * 0.3989422804014327 / sigma;
}

__device__ double
erfc_scaled(double a, double b, double root_two_sigma)
{
  return erfc((a - b) / root_two_sigma);
}

#endif
