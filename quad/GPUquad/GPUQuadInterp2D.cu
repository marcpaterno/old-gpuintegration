#ifndef GPUQUADINTERP2D_H
#define GPUQUADINTERP2D_H

#include "../util/cudaTimerUtil.h"
#include "../util/cudaMemoryUtil.h"

namespace quad{

  

  class Interp2D{

  public:
    Interp2D(){};
    
    cudaTextureObject_t interpT;
    cudaTextureObject_t interpR;
    cudaTextureObject_t interpC;
    int _rows;
    int _cols;
    
    void Allocate(float *arr, int rows, int cols)
    {

      allocText(interpT, arr, rows, cols);
    }
    

    void Initialize(float *xs, float *ys, float *zs, int rows, int cols)
    {
      //rows == 22
      _rows = rows;
      _cols = cols;
      allocText(interpT, zs, rows, cols);
      allocText(interpR, xs, 1, rows);
      allocText(interpC, ys, 1, cols);
      //printf("Irows:%i Icols:%i\n", _rows, _cols);
    }
    
    void Initialize(Interp2D *orig, float *zs)
    {
      allocText(interpT, zs, orig->_rows, orig->_cols);
      _rows = orig->_rows; _cols = orig->_cols;
      interpR = orig->interpR;
      interpC = orig->interpC;
    }
    
    __device__ int SmallNearestN(cudaTextureObject_t array, float target, int size){
      //look for closest value smaller than target or minimum value
      if(tex2D<float>(array,0, 0)> target){	
	return 0;
      }

      //find smallest neighbor, by finding first larger than target, then move one spot back
      for(int i=1; i<size; i++)
      {
	if(tex2D<float>(array, 0, i)>=target){
	  return i-1;
	}
      }
      return size-1; //no appropriate value found, so return the last one
    }


    __device__ float Interp(double y, double x /*, cudaTextureObject_t interpT, cudaTextureObject_t interpC, cudaTextureObject_t interpR, int offset*/){

      float y1, y2, x1, x2;//interpolation points, y-dimension is row dimension
      int y1_index = SmallNearestN(interpR, y, _rows);          //get the index of nearest smaller neighbor

      y1 = tex2D<float>(interpR, 0,   y1_index);
      y2 = tex2D<float>(interpR, 0,   y1_index+1);
      
      int x1_index = SmallNearestN(interpC, x, _cols);
      x1 = tex2D<float>(interpC, 0, x1_index);
      x2 = tex2D<float>(interpC, 0, x1_index+1);
      
      //if(threadIdx.x==0 && blockIdx.x==0)
      //printf("X:%f Y:%f | xIndex:%i yIndex:%i | x1:%f x2:%f y1:%f y2:%f\n",x, y, x1_index, y1_index, x1, x2, y1, y2);

      float q11, q12, q21, q22; //Function values for interpolation points
      int offset = 0;
      q11 = tex2D<float>(interpT, y1_index, x1_index+(offset*INTERP_COLS));
      q12 = tex2D<float>(interpT, y1_index+1, x1_index+(offset*INTERP_COLS));
      q21 = tex2D<float>(interpT ,y1_index, x1_index+1+(offset*INTERP_COLS));
      q22 = tex2D<float>(interpT ,y1_index+1, x1_index+1+(offset*INTERP_COLS));

      double  t1 = (x2-x)/((x2-x1)*(y2-y1));
      double  t2 = (x-x1)/((x2-x1)*(y2-y1));
      return ((q11*(y2-y) + q12*(y-y1))*t1 + (q21*(y2-y)+q22*(y-y1))*t2);
      
    }

    
    __device__ double operator()(double x, double y){

	   return Interp(x, y);      
    }

  };
}

#endif
