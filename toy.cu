#include <array>
#include <experimental/tuple>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <tuple>
#include <type_traits>
#include <utility>

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

__host__ __device__ double
foo(double x, double y, double z)
{

  printf("x:%f y:%f z:%f\n", x, y, z);
  return x + y + z;
}

// class Test {
// public:
// __host__ __device__
// Test()
// {}
// __host__ __device__ // Test(int r)
// {
//   rows = r;
// }
// __host__ __device__ double
// operator()(double x, double y, double z) const
// {
//   printf("rows:%i\n", rows);
//   return 0.0;
// }
// int rows;
// };

namespace gpu {

  namespace detail {
    template <class F, size_t N, std::size_t... I>
    __device__ double
    apply_impl(F&& f,
               gpu::cuArray<double, N> const& data,
               std::index_sequence<I...>)
    {
      return f(data[I]...);
    };
  }

  template <class F, size_t N>
  __device__ double
  // Unsure if we need to pass 'f' by value, for GPU execution
  apply(F&& f, gpu::cuArray<double, N> const& data)
  {
    return detail::apply_impl(
      std::forward<F>(f), data, std::make_index_sequence<N>());
  }
}

__global__ void
toyKernel(gpu::cuArray<double, 3> const* args, double* result)
{
  *result = gpu::apply(foo, *args);
}

double
execute_sample(std::array<double, 3> const& args)
{
  gpu::cuArray<double, 3>* d_args;
  std::size_t const nbytes = args.size() * sizeof(double);
  cudaMalloc((void**)&d_args, nbytes);
  cudaMemcpy(d_args, args.data(), nbytes, cudaMemcpyHostToDevice);
  double* d_answer;
  cudaMalloc((void**)&d_answer, sizeof(double));
  toyKernel<<<1, 1>>>(d_args, d_answer);
  cudaDeviceSynchronize();
  cudaFree(d_args);
  double answer = -1.0;
  cudaMemcpy(&answer, d_answer, sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_answer);
  return answer;
}

int
main()
{
  // Test t;
  std::array<double, 3> cpuarray = {1, 2, 3};
  std::cout << execute_sample(cpuarray) << '\n';
  std::experimental::apply(foo, cpuarray);
  return 0;
}
