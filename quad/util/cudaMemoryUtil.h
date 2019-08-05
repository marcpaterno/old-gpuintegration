#ifndef CUDAMEMORYUTIL_H
#define CUDAMEMORYUTIL_H
#include "cudaDebugUtil.h"
void
allocText(cudaTextureObject_t& dest, float table[], int rows, int cols)
{

  cudaArray* cuArray;
  cudaChannelFormatDesc channelDesc =
    cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  cudaMallocArray(&cuArray, &channelDesc, rows, cols);
  cudaMemcpyToArray(
    cuArray, 0, 0, table, sizeof(float) * rows * cols, cudaMemcpyHostToDevice);

  struct cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = cuArray;

  struct cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.addressMode[0] = cudaAddressModeWrap;
  texDesc.addressMode[1] = cudaAddressModeWrap;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = 0;

  dest = 0;
  cudaCreateTextureObject(&dest, &resDesc, &texDesc, NULL);
}

void
allocText(cudaTextureObject_t& dest,
          float table1[],
          float table2[],
          float table3[],
          float table4[],
          float table5[],
          int rows,
          int cols)
{

  cudaArray* cuArray;
  cudaChannelFormatDesc channelDesc =
    cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  cudaMallocArray(&cuArray, &channelDesc, rows, cols * 5);

  cudaMemcpyToArray(
    cuArray, 0, 0, table1, sizeof(float) * rows * cols, cudaMemcpyHostToDevice);
  cudaMemcpyToArray(cuArray,
                    0,
                    cols,
                    table2,
                    sizeof(float) * rows * cols,
                    cudaMemcpyHostToDevice);
  cudaMemcpyToArray(cuArray,
                    0,
                    cols * 2,
                    table3,
                    sizeof(float) * rows * cols,
                    cudaMemcpyHostToDevice);
  cudaMemcpyToArray(cuArray,
                    0,
                    cols * 3,
                    table4,
                    sizeof(float) * rows * cols,
                    cudaMemcpyHostToDevice);
  cudaMemcpyToArray(cuArray,
                    0,
                    cols * 4,
                    table5,
                    sizeof(float) * rows * cols,
                    cudaMemcpyHostToDevice);

  struct cudaResourceDesc resDesc;
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = cuArray;

  struct cudaTextureDesc texDesc;
  memset(&texDesc, 0, sizeof(texDesc));
  texDesc.addressMode[0] = cudaAddressModeWrap;
  texDesc.addressMode[1] = cudaAddressModeWrap;
  texDesc.filterMode = cudaFilterModePoint;
  texDesc.readMode = cudaReadModeElementType;
  texDesc.normalizedCoords = 0;

  dest = 0;
  cudaCreateTextureObject(&dest, &resDesc, &texDesc, NULL);
}

namespace quad {
  template <typename T>
  class MemoryUtil {};

  template <typename T>
  class HostMemory : public MemoryUtil<T> {

  public:
    void*
    AllocateMemory(void* ptr, size_t n)
    {
      ptr = malloc(n);
      return ptr;
    }

    void
    ReleaseMemory(void* ptr)
    {
      free(ptr);
    }
  };

  template <typename T>
  class DeviceMemory : public MemoryUtil<T> {
  public:
    cudaError_t
    AllocateMemory(void** d_ptr, size_t n)
    {
      return cudaMalloc(d_ptr, n);
    }

    cudaError_t
    ReleaseMemory(void* d_ptr)
    {
      return cudaFree(d_ptr);
    }

    cudaError_t
    SetHeapSize(size_t hSize = (size_t)2 * 1024 * 1024 * 1024)
    {
      return cudaDeviceSetLimit(cudaLimitMallocHeapSize, hSize);
    }

    cudaError_t
    CopyHostToDeviceConstantMemory(const char* d_ptr, void* h_ptr, size_t n)
    {
      return cudaMemcpyToSymbol(d_ptr, h_ptr, n);
    }

    cudaError_t
    CopyHostToDeviceConstantMemory(const void* d_ptr, void* h_ptr, size_t n)
    {
      return cudaMemcpyToSymbol(d_ptr, h_ptr, n);
    }

    //@brief Initialize Device
    cudaError_t
    DeviceInit(int dev = -1, int verbose = 0)
    {
      cudaError_t error = cudaSuccess;

      do {
        int deviceCount;
        error = QuadDebug(cudaGetDeviceCount(&deviceCount));
        if (error)
          break;
        if (deviceCount == 0) {
          fprintf(stderr, "No devices supporting CUDA.\n");
          exit(1);
        }

        if ((dev > deviceCount - 1) || (dev < 0)) {
          dev = 0;
        }

        error = QuadDebug(cudaSetDevice(dev));
        if (error)
          break;

        size_t free_physmem, total_physmem;
        QuadDebugExit(cudaMemGetInfo(&free_physmem, &total_physmem));

        cudaDeviceProp deviceProp;
        error = QuadDebug(cudaGetDeviceProperties(&deviceProp, dev));
        if (error)
          break;

        if (deviceProp.major < 1) {
          fprintf(stderr, "Device does not support CUDA.\n");
          exit(1);
        }

        if (false && verbose) {
          printf("Using device %d: %s (SM%d, %d SMs, %lld free / %lld total MB "
                 "physmem, ECC %s)\n",
                 dev,
                 deviceProp.name,
                 deviceProp.major * 100 + deviceProp.minor * 10,
                 deviceProp.multiProcessorCount,
                 (unsigned long long)free_physmem / 1024 / 1024,
                 (unsigned long long)total_physmem / 1024 / 1024,
                 (deviceProp.ECCEnabled) ? "on" : "off");
          fflush(stdout);
        }

      } while (0);
      return error;
    }
  };

}
#endif
