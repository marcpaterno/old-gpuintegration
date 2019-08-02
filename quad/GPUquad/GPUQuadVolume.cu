struct Volume {

  cudaTextureObject_t high;
  cudaTextureObject_t low;

  __device__ float
  GetHigh(int i)
  {
    return tex2D<float>(high, 0, i);
  }
  __device__ float
  GetLow(int i)
  {
    return tex2D<float>(low, 0, i);
  }

  void
  Initialize(float* l, float* h, int dim)
  {
    allocText(high, h, 1, dim);
    allocText(low, l, 1, dim);
  }
};
