#include "GPUQuadSample.cu"

namespace quad{
	
  template<typename T, typename integT>
    __device__
    void
  INIT_REGION_POOL(T *dRegions, T *dRegionsLength, size_t numRegions, int fEvalPerRegion, const int NDIM, integT *test, Volume *vol = NULL)
 {
   size_t index = blockIdx.x;
		
   if(threadIdx.x == 0)
     {	
       for(int dim = 0; dim < NDIM; ++dim)
	 {	
	   T lower = dRegions[dim * numRegions + index];
	   sRegionPool[threadIdx.x].bounds[dim].lower = 0;
	   sRegionPool[threadIdx.x].bounds[dim].upper = 1;
	   sBound[dim].unScaledLower = lower;
	   sBound[dim].unScaledUpper = lower + dRegionsLength[dim * numRegions + index];
	   sRegionPool[threadIdx.x].div = 0;
	 }
     }
	  
   __syncthreads();
   SampleRegionBlock<T>(0, fEvalPerRegion, NDIM, test, vol);
   __syncthreads();
	  
 }

  template<typename T, typename integT>
    __global__
    void 
  INTEGRATE_GPU_PHASE1(T *dRegions, T *dRegionsLength, size_t numRegions, T *dRegionsIntegral, T *dRegionsError, int *activeRegions, int *subDividingDimension, T epsrel, T epsabs, int fEvalPerRegion, const int NDIM, integT *test, Volume *vol = NULL){
		
   
    T ERR = 0, RESULT = 0;
    int fail = 0;	
    INIT_REGION_POOL(dRegions, dRegionsLength, numRegions, fEvalPerRegion, NDIM, test, vol);
   
       //each block has a region
    if(threadIdx.x == 0){
      ERR = sRegionPool[threadIdx.x].result.err;
      RESULT = sRegionPool[threadIdx.x].result.avg;
      T ratio  = ERR/MaxErr(RESULT, epsrel, epsabs);
      int fourthDiffDim = sRegionPool[threadIdx.x].result.bisectdim;
      dRegionsIntegral[gridDim.x + blockIdx.x] = RESULT;
      dRegionsError[gridDim.x + blockIdx.x] = ERR;
			
      if(ratio > 1){
	fail = 1;
	ERR = 0;
	RESULT = 0;
	//printf("Region Failed\n");
      }
      activeRegions[blockIdx.x] = fail;
      subDividingDimension[blockIdx.x] = fourthDiffDim;
      dRegionsIntegral[blockIdx.x] = RESULT;
      dRegionsError[blockIdx.x]=ERR;    
    }  
  }



  ////PHASE 2 Procedures Starts


	template<typename T>
    __device__ void
    ComputeErrResult(T &ERR, T &RESULT){
    /*sdata[threadIdx.x] = sRegionPool[threadIdx.x].result.err;
    sdata[blockDim.x + threadIdx.x] = sRegionPool[threadIdx.x].result.avg;
    __syncthreads();

    // contiguous range pattern
    for(size_t offset = size / 2; offset > 0; offset >>= 1){
      if(threadIdx.x < offset){
	sdata[threadIdx.x] += sdata[threadIdx.x + offset];
	sdata[blockDim.x + threadIdx.x] += sdata[blockDim.x + threadIdx.x + offset];
      }
      __syncthreads();
    }
    */
		if(threadIdx.x == 0){
			ERR = sRegionPool[threadIdx.x].result.err;
			RESULT = sRegionPool[threadIdx.x].result.avg;
		}
		__syncthreads();
	}
	
  template<typename T, typename integT>
    __device__ int
  INIT_REGION_POOL(T *dRegions, T *dRegionsLength, int *subDividingDimension,  size_t numRegions,  int fEvalPerRegion, const int NDIM, integT *test, Volume *vol = NULL){
				
		size_t intervalIndex = blockIdx.x;
		int idx = 0;
			
		for(; idx < SM_REGION_POOL_SIZE/BLOCK_SIZE; ++idx){
			int index = idx*BLOCK_SIZE + threadIdx.x;
			sRegionPool[index].div = 0;
			sRegionPool[index].result.err = 0;
			sRegionPool[index].result.avg = 0;
			sRegionPool[index].result.bisectdim = 0;
  
			for(int dim = 0; dim < NDIM; ++dim){
				sRegionPool[index].bounds[dim].lower = 0;
				sRegionPool[index].bounds[dim].upper = 0;
			}
		}
		
		int index = idx*BLOCK_SIZE + threadIdx.x;
		
		if(index < SM_REGION_POOL_SIZE){
			sRegionPool[index].div = 0;
			sRegionPool[index].result.err = 0;
			sRegionPool[index].result.avg = 0;
			sRegionPool[index].result.bisectdim = 0;
  
			for(int dim = 0; dim < NDIM; ++dim){
				sRegionPool[index].bounds[dim].lower = 0;
				sRegionPool[index].bounds[dim].upper = 0;
			}      
		}
    
		if(threadIdx.x == 0){
			for(int dim = 0; dim < NDIM; ++dim){
				sRegionPool[threadIdx.x].bounds[dim].lower = 0;
				sRegionPool[threadIdx.x].bounds[dim].upper = 1;
				T lower = dRegions[dim * numRegions + intervalIndex];
				sBound[dim].unScaledLower = lower;
				sBound[dim].unScaledUpper = lower + dRegionsLength[dim * numRegions + intervalIndex];
			}    
		}
    
		__syncthreads();
		SampleRegionBlock<T>(0, fEvalPerRegion, NDIM, test, vol);
    
		if(threadIdx.x == 0){
			gPool = (Region *)malloc(sizeof(Region)*(SM_REGION_POOL_SIZE/2));
			gRegionPoolSize = (SM_REGION_POOL_SIZE/2);//BLOCK_SIZE;
		}

		__syncthreads();
    
		for(idx = 0; idx < (SM_REGION_POOL_SIZE/2)/BLOCK_SIZE; ++idx){
			int index = idx*BLOCK_SIZE + threadIdx.x;
			gRegionPos[index]=index;
			gPool[index]=sRegionPool[index];
		}
		index = idx*BLOCK_SIZE + threadIdx.x;
		if(index < (SM_REGION_POOL_SIZE/2)){
			gRegionPos[index]=index;
			gPool[index]=sRegionPool[index];
		}
		return 1;
	}
  
	template <class T>
    __device__
    void swap ( T& a, T& b ){
    T c(a);
    a=b; b=c;
  }
  
	template<typename T>
    __device__
    void
    INSERT_GLOBAL_STORE(Region *sRegionPool, Region *gRegionPool){
		if(threadIdx.x == 0){
			gPool = (Region *)malloc(sizeof(Region)*(gRegionPoolSize+(SM_REGION_POOL_SIZE/2)));
		}
		__syncthreads();
    
 
		//Copy existing global regions into newly allocated spaced
		int iterationsPerThread = 0;
		for(iterationsPerThread = 0; iterationsPerThread < gRegionPoolSize/BLOCK_SIZE; ++iterationsPerThread){
			size_t dataIndex = iterationsPerThread*BLOCK_SIZE+threadIdx.x;
			gPool[dataIndex] = gRegionPool[dataIndex];
		}
    
		size_t dataIndex = iterationsPerThread*BLOCK_SIZE+threadIdx.x;
		if(dataIndex < gRegionPoolSize){
			gPool[dataIndex] = gRegionPool[dataIndex];
		}

		//Fill the previous occupied postion in global memory by half of shared memory regions
		for(iterationsPerThread = 0; iterationsPerThread < (SM_REGION_POOL_SIZE/2)/BLOCK_SIZE; ++iterationsPerThread){
			int index = iterationsPerThread*BLOCK_SIZE+threadIdx.x;
			gPool[gRegionPos[index]] = sRegionPool[index];
			gPool[gRegionPoolSize + index] = sRegionPool[(SM_REGION_POOL_SIZE/2) + index];

		}
		int index = iterationsPerThread*BLOCK_SIZE+threadIdx.x;
		if(index < (SM_REGION_POOL_SIZE/2)){
			int index = iterationsPerThread*BLOCK_SIZE+threadIdx.x;
			gPool[gRegionPos[index]] = sRegionPool[index];
			gPool[gRegionPoolSize + index] = sRegionPool[(SM_REGION_POOL_SIZE/2) + index];
		}

		//gPool[gRegionPos[threadIdx.x]] = sRegionPool[threadIdx.x];
		//gPool[gRegionPoolSize + threadIdx.x] = sRegionPool[BLOCK_SIZE + threadIdx.x];

		__syncthreads();
		if(threadIdx.x == 0){
			gRegionPoolSize = gRegionPoolSize+(SM_REGION_POOL_SIZE/2);
			free(gRegionPool);
		}
		__syncthreads();
		gRegionPool = gPool;
		// gSize += BLOCK_SIZE;

		//return gSize;
  }

	template<typename T>
    __device__ void
    EXTRACT_MAX(T *serror, size_t *serrorPos, size_t gSize)
	{
		for(size_t offset = gSize/2; offset > 0; offset >>= 1 ){
			int idx = 0;
			for(idx = 0; idx < offset/BLOCK_SIZE; ++idx){
				size_t index = idx*BLOCK_SIZE+threadIdx.x;
				if(index < offset){
					if(serror[index] < serror[index+offset]){
						swap(serror[index], serror[index+offset]);
						swap(serrorPos[index], serrorPos[index+offset]);
					}
					//printf("%ld %ld\n",index, index+offset);
				}
			}
			size_t index = idx*BLOCK_SIZE+threadIdx.x;
			if(index < offset){
				if(serror[index] < serror[index+offset]){
					swap(serror[index], serror[index+offset]);
					swap(serrorPos[index], serrorPos[index+offset]);
				}
			//printf("%ld %ld\n",index, index+offset);
			}
			__syncthreads();
		}
	}
	
	template<typename T>
    __device__ void
    EXTRACT_TOPK(Region *sRegionPool, Region *gRegionPool){

    //Comment 3 instructions these section if you are directly using new shared memory instead of reusing shared memory
    
    T *sarray = (T *)&sRegionPool[0];

    if(threadIdx.x == 0){
      //T *sarray = (T *)&sRegionPool[0];

      if((gRegionPoolSize*sizeof(T) + gRegionPoolSize*sizeof(size_t)) < sizeof(Region) * SM_REGION_POOL_SIZE){
        serror = &sarray[0];
	//TODO:Size of sRegionPool vs sarray constrain
	serrorPos = (size_t *)&sarray[gRegionPoolSize];
      }else{
	serror = (T *)malloc(sizeof(T) * gRegionPoolSize);
	serrorPos = (size_t *)malloc(sizeof(size_t) * gRegionPoolSize);
      }
    }
    __syncthreads();

    int offset = 0;
    for(offset = 0; (offset < MAX_GLOBALPOOL_SIZE/BLOCK_SIZE) && (offset < gRegionPoolSize/BLOCK_SIZE); offset++){
      size_t regionIndex = offset*BLOCK_SIZE + threadIdx.x;
      serror[regionIndex] = gRegionPool[regionIndex].result.err;
      serrorPos[regionIndex] = regionIndex;
    }
    size_t regionIndex = offset*BLOCK_SIZE + threadIdx.x;
    if(regionIndex < gRegionPoolSize){
      serror[regionIndex] = gRegionPool[regionIndex].result.err;
      serrorPos[regionIndex] = regionIndex;
    }

    __syncthreads();
    for(int k = 0; k < (SM_REGION_POOL_SIZE/2); ++k){
      EXTRACT_MAX<T>(&serror[k], &serrorPos[k], gRegionPoolSize-k);
    }
    
    int iterationsPerThread = 0;
    for(iterationsPerThread = 0; iterationsPerThread < (SM_REGION_POOL_SIZE/2)/BLOCK_SIZE; ++iterationsPerThread){
      int index = iterationsPerThread*BLOCK_SIZE + threadIdx.x;
      size_t pos = serrorPos[index];
      gRegionPos[index] = pos;
    }
    int index = iterationsPerThread*BLOCK_SIZE + threadIdx.x;
    if(index < (SM_REGION_POOL_SIZE/2)){
      size_t pos = serrorPos[index];
      gRegionPos[index] = pos;
    }

    //Old
    //size_t pos = serrorPos[threadIdx.x];
    //gRegionPos[threadIdx.x] = pos;
    __syncthreads();

    if(threadIdx.x == 0){
      //sRegionPool = (Region *)sarray;
      if(2*gRegionPoolSize*sizeof(T) >= sizeof(Region) * SM_REGION_POOL_SIZE){
	free(serror);
	free(serrorPos);
      }
    }
    __syncthreads();

    /*if((2*gRegionPoolSize*sizeof(T) >= sizeof(Region) * SM_REGION_POOL_SIZE) && threadIdx.x == 0){
      free(serror);
      free(serrorPos);
      }*/
    
    //Copy top K into SM and reset the remaining
    for(iterationsPerThread = 0; iterationsPerThread < (SM_REGION_POOL_SIZE/2)/BLOCK_SIZE; ++iterationsPerThread){
      int index = iterationsPerThread*BLOCK_SIZE + threadIdx.x;
      sRegionPool[index] = gPool[gRegionPos[index]];
      sRegionPool[(SM_REGION_POOL_SIZE/2) + index].result.err = -INFTY;
      sRegionPool[(SM_REGION_POOL_SIZE/2) + index].result.avg = 0;
      sRegionPool[(SM_REGION_POOL_SIZE/2) + index].div = 0;
    }
    index = iterationsPerThread*BLOCK_SIZE + threadIdx.x;
    if(index < (SM_REGION_POOL_SIZE/2)){
      sRegionPool[index] = gPool[gRegionPos[index]];
      sRegionPool[(SM_REGION_POOL_SIZE/2) + index].result.err = -INFTY;
      sRegionPool[(SM_REGION_POOL_SIZE/2) + index].result.avg = 0;
      sRegionPool[(SM_REGION_POOL_SIZE/2) + index].div = 0;
    }

    //Old
    //sRegionPool[threadIdx.x] = gPool[pos];

    //sRegionPool[BLOCK_SIZE+threadIdx.x].result.err = -INFTY;
    //sRegionPool[BLOCK_SIZE+threadIdx.x].result.avg = 0;
    //sRegionPool[BLOCK_SIZE+threadIdx.x].div = 0;


  }


  template<typename T>
    __device__ size_t
    EXTRACT_MAX(Region *sRegionPool, Region *gRegionPool, size_t sSize)
	{
		//If SharedPool is full
		if(sSize == SM_REGION_POOL_SIZE){
			INSERT_GLOBAL_STORE<T>(sRegionPool, gRegionPool);
			gRegionPool = gPool;
			EXTRACT_TOPK<T>(sRegionPool, gRegionPool);
			sSize = (SM_REGION_POOL_SIZE/2);
			__syncthreads();
		//    return sSize;
		}
    
		for(size_t offset = (SM_REGION_POOL_SIZE/2); offset > 0; offset >>= 1 ){
			int idx = 0;
			for(idx = 0; idx < offset/BLOCK_SIZE; ++idx){
				size_t index = idx*BLOCK_SIZE+threadIdx.x;	
				if(index < offset){
					Region *r1 = &sRegionPool[index];
					Region *r2 = &sRegionPool[index + offset];
					if(r1->result.err < r2->result.err ){
						swap<Region>(sRegionPool[index],sRegionPool[offset + index]);
					}
				}
			}
			size_t index = idx*BLOCK_SIZE+threadIdx.x;
			if(index < offset){
				Region *r1 = &sRegionPool[index];
				Region *r2 = &sRegionPool[index + offset];
				if(r1->result.err < r2->result.err ){
					swap<Region>(sRegionPool[index],sRegionPool[offset + index]);
				}
			}
			__syncthreads();
    }	

    return sSize;
  }


  template<typename T, typename integT>
    __global__ void
  BLOCK_INTEGRATE_GPU_PHASE2(T *dRegions, T *dRegionsLength, size_t numRegions, T *dRegionsIntegral, T *dRegionsError, int *dRegionsNumRegion, int *activeRegions, int *subDividingDimension, T epsrel, T epsabs, int fEvalPerRegion , const int NDIM, integT *test, Volume *vol = NULL){
    
    
    Region *gRegionPool = 0;
    int sRegionPoolSize = INIT_REGION_POOL<T>(dRegions, dRegionsLength, subDividingDimension,  numRegions, fEvalPerRegion, NDIM, test, vol);
    	
    ComputeErrResult<T>(ERR, RESULT);
  
    //TODO : May be redundance sync
    __syncthreads();
    int nregions = sRegionPoolSize;
		
    for(; (nregions <= MAX_GLOBALPOOL_SIZE)&&(nregions == 1 || ERR > MaxErr(RESULT, epsrel, epsabs)); ++nregions ){
      gRegionPool = gPool;
      sRegionPoolSize = EXTRACT_MAX<T>(sRegionPool, gRegionPool, sRegionPoolSize);
			
      Region *RegionLeft, *RegionRight;
      Result result;
      if(threadIdx.x == 0){
	Bounds *bL, *bR;
	Region *R = &sRegionPool[0];
	result.err = R->result.err;
	result.avg = R->result.avg;
	result.bisectdim = R->result.bisectdim;

	int bisectdim = result.bisectdim;

	RegionLeft = R;
	RegionRight = &sRegionPool[sRegionPoolSize];

	bL = &RegionLeft->bounds[bisectdim];
	bR = &RegionRight->bounds[bisectdim];
				
	RegionRight->div = ++RegionLeft->div;
	for(int dim = 0; dim < NDIM; ++dim){
	  RegionRight->bounds[dim].lower = RegionLeft->bounds[dim].lower;
	  RegionRight->bounds[dim].upper = RegionLeft->bounds[dim].upper;
	}
	//Subdivide the chosen axis
	bL->upper = bR->lower = 0.5 * (bL->lower + bL->upper);	
      }
      sRegionPoolSize++;
      __syncthreads();
      SampleRegionBlock<T>(0,fEvalPerRegion, NDIM, test, vol);

      __syncthreads();
      SampleRegionBlock<T>(sRegionPoolSize-1, fEvalPerRegion, NDIM, test, vol);

      __syncthreads();
      
      if(threadIdx.x == 0){
	Result *rL = &RegionLeft->result;
	Result *rR = &RegionRight->result;

	T diff = rL->avg + rR->avg - result.avg;
	diff = fabs(.25*diff);
	T err = rL->err + rR->err;
	if( err > 0 ) {
	  T c = 1 + 2*diff/err;
	  rL->err *= c;
	  rR->err *= c;
	}
	rL->err += diff;
	rR->err += diff;

	ERR += rL->err + rR->err - result.err;
	RESULT +=  rL->avg + rR->avg - result.avg;

      }
      __syncthreads();
    }
    
    if(threadIdx.x == 0){
      int isActive = ERR > MaxErr(RESULT, epsrel, epsabs);
      if(/*(nregions > MAX_GLOBALPOOL_SIZE) || isActive || */ERR > (1e+10)){
	//printf("%ld %.16lf %.16lf %ld\n",(size_t)blockIdx.x, RESULT, ERR, (size_t)nregions);
	RESULT = 0.0;
	ERR = 0.0;
	isActive = 1;
      }
      //printf("%ld %.16lf %.16lf %ld\n",(size_t)blockIdx.x, RESULT, ERR, (size_t)nregions);
      activeRegions[blockIdx.x] = isActive;
      dRegionsIntegral[blockIdx.x] = RESULT;
      dRegionsError[blockIdx.x] = ERR;
      dRegionsNumRegion[blockIdx.x] = nregions;
      free(gPool);
    }
  }
    

}
