#include <array>
#include <experimental/tuple>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <tuple>
#include <type_traits>

using namespace std;

namespace gpu {

  template <typename T, size_t s>
  class cuArray {

  public:
    __host__ __device__ const T*
    begin() const
    {
      return &data[0];
    }
    __host__ __device__ const T*
    end() const
    {
      return (&data[0] + s);
    }
    __host__ __device__ constexpr size_t
    size() const
    {
      return s;
    }
    __host__ __device__ T& operator[](const size_t i) { return data[i]; }
    __host__ __device__ T const& operator[](size_t i) const { return data[i]; }
    T data[s];
  };
};

namespace std {
  template <size_t index, typename T, size_t s>
  /*__device__ __host__*/
  constexpr T&
  get(gpu::cuArray<T, s> const& x) noexcept
  {
    return x[index];
  }

  template <class T, size_t N>
  class tuple_size<gpu::cuArray<T, N>>
    : public std::integral_constant<size_t, N> {};

};

__host__ __device__ double
foo(double x, double y, double z)
{
  printf("Second Parameter %f\n, y");
  return 0;
}

__global__ void
kernel()
{}

int
main()
{
  // test cuArray
  gpu::cuArray<double, 3> carr;
  carr[0] = 1;
  carr[1] = 2;
  carr[2] = 3;
  std::experimental::apply(foo, carr);

  // kernel<<<1,1>>>();
  cudaDeviceSynchronize();
  return 0;
}
