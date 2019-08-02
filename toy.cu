#include <array>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <tuple>
//#include <experimental/tuple>
#include <type_traits>
#include <utility>
// using namespace std;
#include <iostream>

namespace gpu {
  template <typename T, std::size_t s>
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
    __host__ __device__ constexpr std::size_t
    size() const
    {
      return s;
    }
    __host__ __device__ T& operator[](std::size_t i) { return data[i]; }
    __host__ __device__ T const& operator[](std::size_t i) const
    {
      return data[i];
    }
    T data[s];
  };
};

__host__ __device__ void
foo(double x, double y, double z)
{
  printf("x:%f y:%f z:%f\n", x, y, z);
}

namespace std {
  template <std::size_t index, typename T, std::size_t s>
  __device__ __host__ constexpr T&
  get(gpu::cuArray<T, s>& x) noexcept
  {
    return x[index];
  }

  template <std::size_t index, typename T, std::size_t s>
  __device__ __host__ constexpr T const&
  get(gpu::cuArray<T, s> const& x) noexcept
  {
    return x[index];
  }

  template <std::size_t index, typename T, std::size_t s>
  __device__ __host__ constexpr T&&
  get(gpu::cuArray<T, s>&& x) noexcept
  {
    return std::move(x[index]);
  }

  template <class T, std::size_t N>
  class tuple_size<gpu::cuArray<T, N>>
    : public std::integral_constant<std::size_t, N> {};
};
#include <experimental/tuple>

class Test {
public:
  __host__ __device__
  Test()
  {}
  __host__ __device__
  Test(int r)
  {
    rows = r;
  }
  __host__ __device__ double
  operator()(double x, double y, double z) const
  {
    printf("rows:%i\n", rows);
  }
  int rows;
};

__global__ void
kernelApply()
{
  gpu::cuArray<double, 3> arr = {1, 2, 3};
  std::apply(foo, arr);
}

int
main()
{
  Test t;
  gpu::cuArray<double, 3> arr = {1, 2, 3};
  std::array<double, 3> arr2 = {1, 2, 3};
  std::apply(foo, arr);
  cudaDeviceSynchronize();
  return 0;
}
