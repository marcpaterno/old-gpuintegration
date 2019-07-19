#include "GPUQuadPhases.cu"
#include <iostream> 
#include "GPUQuadInterp2D.cu"


__global__ void kernel(cudaTextureObject_t texObj, cudaTextureObject_t texObj2)
	{
		printf("hello[0][0]:%f\n", 			tex2D<float>(texObj, 0, 0));
		printf("hello[0][1]:%f\n", 			tex2D<float>(texObj, 0, 1));
		printf("hello[1][0]:%f\n", 			tex2D<float>(texObj, 0, 2));
		printf("hello[1][1]:%f\n", 			tex2D<float>(texObj, 0, 3));	      
		printf("hey[0][0]:%f\n", 			tex2D<float>(texObj2, 0, 0));
	}

namespace quad{
  __constant__ size_t dFEvalPerRegion;
  template<typename T> 
    __global__ void generateInitialRegions(T *dRegions, T *dRegionsLength, size_t numRegions,  T *newRegions, 
										   T *newRegionsLength, size_t newNumOfRegions, int numOfDivisionsPerRegionPerDimension, int NDIM)
	{
	  extern __shared__ T slength[];
	  size_t threadId = blockIdx.x * blockDim.x + threadIdx.x;
	
	  if(threadIdx.x < NDIM){
	    slength[threadIdx.x] = dRegionsLength[threadIdx.x]/numOfDivisionsPerRegionPerDimension;
	    //printf("thread:$i slength:%f\n",  slength[threadIdx.x]);
	  }
	  __syncthreads();
		
	  if(threadId < newNumOfRegions){
	    size_t interval_index = threadId / pow((T)numOfDivisionsPerRegionPerDimension, (T)NDIM);
	    size_t local_id = threadId %  (size_t)pow((T)numOfDivisionsPerRegionPerDimension, (T)NDIM);
	    for(int dim = 0; dim < NDIM; ++dim){
	      size_t id = (size_t)(local_id/pow((T)numOfDivisionsPerRegionPerDimension, (T)dim)) % numOfDivisionsPerRegionPerDimension;
	      newRegions[newNumOfRegions*dim + threadId] = dRegions[numRegions*dim + interval_index] + id*slength[dim];
	      newRegionsLength[newNumOfRegions*dim + threadId] = slength[dim];
	      //printf("Thread:%lu newRegions:(dim:%i)%f length:%f\n", threadId , dim, newRegions[newNumOfRegions*dim + threadId], slength[dim]);
	    }
	  }
	}
	
  template<typename T> 
    __global__
    void alignRegions(T *dRegions, T *dRegionsLength, int *activeRegions, int *subDividingDimension, 
						int *scannedArray, T *newActiveRegions, T *newActiveRegionsLength, 
		      int *newActiveRegionsBisectDim, size_t numRegions, size_t newNumRegions, int numOfDivisionOnDimension, int NDIM)
	{
		size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
		if(tid <  numRegions && activeRegions[tid] == 1){
			
			size_t interval_index = scannedArray[tid];
			
			for(int i = 0 ; i < NDIM; ++i){
				newActiveRegions[i*newNumRegions + interval_index] = dRegions[i*numRegions + tid];
				newActiveRegionsLength[i*newNumRegions + interval_index] = dRegionsLength[i*numRegions + tid];
			}
			
			for(int i = 0; i < numOfDivisionOnDimension; ++i){
				newActiveRegionsBisectDim[i*newNumRegions + interval_index] = subDividingDimension[tid];
			}
		}
	}

	template<typename T>
    __global__ void
	divideIntervalsGPU(T *genRegions, T *genRegionsLength, T *activeRegions, T *activeRegionsLength, int *activeRegionsBisectDim, size_t numActiveRegions, int numOfDivisionOnDimension, int NDIM)
	{
		size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
		if(tid < numActiveRegions)
		{
			int bisectdim = activeRegionsBisectDim[tid];
			size_t data_size = numActiveRegions*numOfDivisionOnDimension;
			
			//if(threadIdx.x == 0 && blockIdx.x==0)
				//printf("Bisect Dim:%i # of Divisions:%i\n", bisectdim, numOfDivisionOnDimension);
			
			for(int i = 0; i < numOfDivisionOnDimension; ++i){
				for(int dim = 0; dim < NDIM; ++dim){
					genRegions[i*numActiveRegions + dim*data_size + tid]=activeRegions[dim*numActiveRegions + tid];
					genRegionsLength[i*numActiveRegions + dim*data_size + tid]=activeRegionsLength[dim*numActiveRegions + tid];
					
					//if(threadIdx.x == 0 && blockIdx.x==0)
						//printf("genRegions:%f Length:%f\n", activeRegions[dim*numActiveRegions + tid], activeRegionsLength[dim*numActiveRegions + tid]);
				}
			}
			
			//if(threadIdx.x == 0 && blockIdx.x==0)
				//printf("=====================\n");
			for(int i = 0; i < numOfDivisionOnDimension; ++i){
				T interval_length = activeRegionsLength[bisectdim*numActiveRegions + tid]/numOfDivisionOnDimension;
				genRegions[bisectdim*data_size + i*numActiveRegions + tid] = activeRegions[bisectdim*numActiveRegions + tid] + i*interval_length;
				genRegionsLength[i*numActiveRegions + bisectdim*data_size + tid] = interval_length;
				
				//if(threadIdx.x == 0 && blockIdx.x==0)
						//printf("genRegions:%f Length:%f\n", activeRegions[bisectdim*numActiveRegions + tid] + i*interval_length, interval_length);
			}
		}	
	}
	
	template<typename T>
    class GPUKernelCuhre{
		
	  T *dRegions;
	  T *dRegionsLength;
	  T *hRegions;
	  T *hRegionsLength;
		
	  cudaTextureObject_t interpT;
	  cudaTextureObject_t interpR;
	  cudaTextureObject_t interpC;
		
	  int NDIM, KEY, VERBOSE;
	  size_t numRegions, numFunctionEvaluations;
	  size_t fEvalPerRegion;
	  HostMemory<T> Host;
	  DeviceMemory<T> Device;
	  QuadRule<T> Rule;
	  int NUM_DEVICES ;
	  //Debug Msg
	  char msg[256];
		
	  std::ostream &log;
	public:
		
	  GPUKernelCuhre(std::ostream &logerr=std::cout):log(logerr){
	    numRegions = 0;
	    numFunctionEvaluations = 0;
	    NDIM = 0;
	    KEY = 0;
	  }
		
	  ~GPUKernelCuhre(){
	    if(VERBOSE){
	      sprintf(msg, "GPUKerneCuhre Destructur");
	      Println(log, msg);
	    }
		
	    QuadDebug(Device.ReleaseMemory(dRegions));
	    QuadDebug(Device.ReleaseMemory(dRegionsLength));
	    Host.ReleaseMemory(hRegions);
	    Host.ReleaseMemory(hRegionsLength);
	    QuadDebug(cudaDeviceReset());
	    if(DIM > 8)
	      QuadDebug(Device.ReleaseMemory(gpuGenPos));
	  }
    
	  size_t getNumActiveRegions(){
	    return numRegions;
	  }	

	  void setRegionsData(T *data, size_t size){
	    hRegions = &data[0];
	    hRegionsLength = &data[size*NDIM];
	    numRegions = size;
	  }

	  T *getRegions(size_t size, int startIndex){
	    T *newhRegionsAndLength = 0;
	    newhRegionsAndLength = (T *)Host.AllocateMemory(&newhRegionsAndLength, 2*sizeof(T)*size*NDIM);
	    T *newhRegions = &newhRegionsAndLength[0], *newhRegionsLength = &newhRegionsAndLength[size*NDIM];
	    //NOTE:Copy order is important
	    for(int dim = 0; dim < NDIM; ++dim){
	      QuadDebug(cudaMemcpy(newhRegions + dim * size, dRegions + dim * numRegions + startIndex, sizeof(T) * size, cudaMemcpyDeviceToHost));
	      QuadDebug(cudaMemcpy(newhRegionsLength + dim * size, dRegionsLength + dim * numRegions + startIndex, sizeof(T) * size, cudaMemcpyDeviceToHost));
	    } 
	    return newhRegionsAndLength;
	  }

	  void InitGPUKernelCuhre(const int dim, int key, int verbose, int numDevices = 1){
	    QuadDebug(cudaDeviceReset());
	    NDIM = dim;
	    KEY = key;
	    
	    VERBOSE = verbose;
	    NUM_DEVICES = numDevices;
	    fEvalPerRegion = (1 + 2*NDIM + 2*NDIM + 2*NDIM + 2*NDIM + 2*NDIM*(NDIM - 1) + 4*NDIM*(NDIM - 1) + 4*NDIM*(NDIM - 1)*(NDIM - 2)/3 + (1 << NDIM));
	                       
	    QuadDebug(cudaMemcpyToSymbol(dFEvalPerRegion, &fEvalPerRegion, sizeof(size_t), 0, cudaMemcpyHostToDevice));
	    Rule.Init(NDIM, fEvalPerRegion, KEY, VERBOSE);
	    QuadDebug(Device.SetHeapSize());
	    // GenInterpolation();
	  }
		
	  //@brief Template function to display GPU device array variables
	  template <class K> 
	  void display(K *array, size_t size){
	    K *tmp = (K *)malloc(sizeof(K) * size);
	    cudaMemcpy(tmp, array, sizeof(K)*size, cudaMemcpyDeviceToHost);
	    for(int i = 0 ; i < size; ++i){
	      printf("%.20lf \n",(T)tmp[i]);
	    }
	  }

	  void GenerateInitialRegions(){
			
	    hRegions = (T *)Host.AllocateMemory(&hRegions, sizeof(T)*NDIM); //as many regions as dimensions for now
	    hRegionsLength = (T *)Host.AllocateMemory(&hRegionsLength, sizeof(T)*NDIM); 
		  
	    for(int dim = 0 ; dim < NDIM; ++dim)
	      {
		hRegions[dim] = 0;
#if GENZ_TEST == 1
		hRegionsLength[dim] = b[dim];
#else
		hRegionsLength[dim] = 1;	
#endif
	      }
		  
	    QuadDebug(Device.AllocateMemory((void**)&dRegions, sizeof(T)*NDIM));
	    QuadDebug(Device.AllocateMemory((void**)&dRegionsLength, sizeof(T)*NDIM));

	    QuadDebug(cudaMemcpy(dRegions, hRegions, sizeof(T) * NDIM, cudaMemcpyHostToDevice));
	    QuadDebug(cudaMemcpy(dRegionsLength, hRegionsLength, sizeof(T) * NDIM, cudaMemcpyHostToDevice));

	    size_t numThreads = 512;
	    size_t numOfDivisionPerRegionPerDimension = 4;
	    if(NDIM == 5 )
	      numOfDivisionPerRegionPerDimension = 2;
	    if(NDIM == 6 )
	      numOfDivisionPerRegionPerDimension = 2;
	    if(NDIM == 7 )
	      numOfDivisionPerRegionPerDimension = 2;
	    if(NDIM > 7 )
	      numOfDivisionPerRegionPerDimension = 2;
	    if(NDIM > 10 )
	      numOfDivisionPerRegionPerDimension = 1;
		  
	    size_t numBlocks = (size_t)ceil(pow((T)numOfDivisionPerRegionPerDimension, (T)NDIM) / numThreads); //multiple subregions per block
	    numRegions = (size_t)pow((T)numOfDivisionPerRegionPerDimension, (T)NDIM);
		  
	    T *newRegions = 0;
	    T *newRegionsLength = 0;
		  
	    QuadDebug(Device.AllocateMemory((void **)&newRegions, sizeof(T)*numRegions*NDIM));
	    QuadDebug(Device.AllocateMemory((void **)&newRegionsLength, sizeof(T)*numRegions*NDIM));
		  
	    //512 threads per block, we want 256 regions -> 1 block for now
	    generateInitialRegions<T><<<numBlocks, numThreads, NDIM*sizeof(T)>>>(dRegions, dRegionsLength, 1, newRegions, newRegionsLength, 
										 numRegions, numOfDivisionPerRegionPerDimension, NDIM);
	    QuadDebug(Device.ReleaseMemory((void *)dRegions));  //with the new regions generated, we dont' need the  memory of the original ones
	    QuadDebug(Device.ReleaseMemory((void *)dRegionsLength));
		  
	    dRegions = newRegions; //but we still like the variable names
	    dRegionsLength = newRegionsLength;
	    QuadDebug(cudaMemcpy(dRegions, newRegions, sizeof(T) * numRegions * NDIM, cudaMemcpyDeviceToDevice));
	    QuadDebug(cudaMemcpy(dRegionsLength, newRegionsLength, sizeof(T) * numRegions * NDIM, cudaMemcpyDeviceToDevice));
			
	    //QuadDebug(Device.MemoryFree((void **)newRegions));
	    //QuadDebug(Device.MemoryFree((void **)newRegionsLength));
	  }
		 
	  void GenerateActiveIntervals(int *activeRegions, int *subDividingDimension){
			
	    int *scannedArray = 0;
	    QuadDebug(Device.AllocateMemory((void **)&scannedArray, sizeof(int)*numRegions));

	    thrust::device_ptr<int> d_ptr 		= thrust::device_pointer_cast(activeRegions);
	    thrust::device_ptr<int> scan_ptr 	= thrust::device_pointer_cast(scannedArray);
	    thrust::exclusive_scan(d_ptr, d_ptr + numRegions, scan_ptr); //creates running sum, with scan_ptr[0]=0, for the scannedArray
			
	    int last_element;
	    size_t numActiveRegions = 0;
			
	    QuadDebug(cudaMemcpy(&last_element, 	activeRegions 	+ numRegions - 1, 		sizeof(int), 	cudaMemcpyDeviceToHost));
	    QuadDebug(cudaMemcpy(&numActiveRegions,	scannedArray 	+ numRegions - 1, 		sizeof(int), 	cudaMemcpyDeviceToHost));
		  
	    if(last_element == 1)
	      numActiveRegions++;
		  
	    //std::cout << "GenerateActiveIntervals: 	numActiveRegions:" << numActiveRegions << "\n";

	    if(numActiveRegions > 0){

	      int numOfDivisionOnDimension = 2;
				
	      /*if(numActiveRegions < (1 << 7)){
		numOfDivisionOnDimension = (1<<11)/numActiveRegions;
		numOfDivisionOnDimension = 1 << ((int)ceil(log(numOfDivisionOnDimension)/log(2))-1);
		}else if(numActiveRegions < (1<<10)){
		numOfDivisionOnDimension = 8;
		}else if(numActiveRegions < (1<<12)){
		numOfDivisionOnDimension = 4;
		}else{
		numOfDivisionOnDimension = 2;
		}*/
		
	      if(VERBOSE){
		sprintf(msg, "\nComputing NumOfDivisionsOnDimension\n\t#. of Active Regions\t\t: %ld\n\tDivision on dimension\t\t: %ld division", numActiveRegions, numOfDivisionOnDimension);
		Println(log, msg);
	      }
		
	      int *newActiveRegionsBisectDim = 0;
	      T *newActiveRegions = 0, *newActiveRegionsLength = 0;
				
	      cudaMalloc((void **)&newActiveRegions, sizeof(T) *  numActiveRegions * NDIM );
	      cudaMalloc((void **)&newActiveRegionsLength, sizeof(T) *  numActiveRegions * NDIM);
	      cudaMalloc((void **)&newActiveRegionsBisectDim, sizeof(int) * numActiveRegions * numOfDivisionOnDimension);

	      size_t numThreads = BLOCK_SIZE;
	      size_t numBlocks = numRegions/numThreads + ((numRegions%numThreads)?1:0);
		
	      if(VERBOSE){
		Println(log, "\nCalling GPU Function align_intervals");
		sprintf(msg, "\n\t# of input intervals\t\t: %ld\n\t#. of Active Intervals\t\t: %ld\n\t#. of Thread Blocks\t\t: %ld\n\t#. of Threads per Blocks\t: %ld\n",numRegions, numActiveRegions, numBlocks, numThreads);
		Println(log, msg);
	      }

	      alignRegions<T><<<numBlocks, numThreads>>>(dRegions, dRegionsLength, activeRegions, subDividingDimension, scannedArray, newActiveRegions, newActiveRegionsLength, newActiveRegionsBisectDim, numRegions, numActiveRegions, numOfDivisionOnDimension, NDIM);
		
	      if(VERBOSE){
		Println(log, "\nCalling GPU Function divideIntervalsGPU");
		sprintf(msg, "\n\t# of input intervals\t\t: %ld\n\t#. of division on dimension\t: %ld\n\t#. of Thread Blocks\t\t: %ld\n\t#. of Threads per Blocks\t: %ld",numActiveRegions, numOfDivisionOnDimension, numBlocks, numThreads);
		Println(log, msg);
	      }
			
	      T *genRegions = 0, *genRegionsLength = 0;
	      numBlocks = numActiveRegions/numThreads + ((numActiveRegions%numThreads)?1:0);
	      QuadDebug(cudaMalloc((void **)&genRegions, sizeof(T) * numActiveRegions * NDIM * numOfDivisionOnDimension));
	      QuadDebug(cudaMalloc((void **)&genRegionsLength, sizeof(T) * numActiveRegions * NDIM * numOfDivisionOnDimension));
				
	      divideIntervalsGPU<T><<<numBlocks, numThreads>>>(genRegions, genRegionsLength, newActiveRegions, newActiveRegionsLength, newActiveRegionsBisectDim, numActiveRegions, numOfDivisionOnDimension, NDIM);
				
	      QuadDebug(Device.ReleaseMemory(newActiveRegions));
	      QuadDebug(Device.ReleaseMemory(newActiveRegionsLength));
	      QuadDebug(Device.ReleaseMemory(newActiveRegionsBisectDim));
			
	      numRegions = numActiveRegions * numOfDivisionOnDimension;
	      QuadDebug(Device.ReleaseMemory((void *)dRegions));
	      QuadDebug(Device.ReleaseMemory((void *)dRegionsLength));
	      QuadDebug(Device.ReleaseMemory((void *)scannedArray));
		
	      dRegions = genRegions;
	      dRegionsLength = genRegionsLength;
	      //TODO: throws error
	      //QuadDebug(cudaMemcpy(dRegions, genRegions, sizeof(T) * numRegions * NDIM, cudaMemcpyDeviceToDevice));
	      //QuadDebug(cudaMemcpy(dRegionsLength, genRegionsLength, sizeof(T) * numRegions * NDIM, cudaMemcpyDeviceToDevice));
	    }
	    else{
	      numRegions = 0;
	    }
		  
	  }
		
	  template<typename integT>
	  void FirstPhaseIteration(T epsrel, T epsabs, T &integral, T &error, size_t &nregions, size_t &neval, integT *dtest, Volume *vol = NULL){
	    
	    size_t numThreads = BLOCK_SIZE;
	    size_t numBlocks = numRegions;
	    T *dRegionsError = 0, *dRegionsIntegral = 0;
	    QuadDebug(Device.AllocateMemory((void **)&dRegionsIntegral, sizeof(T)*numRegions*2));
	    QuadDebug(Device.AllocateMemory((void **)&dRegionsError, sizeof(T)*numRegions*2));
  
	    int *activeRegions = 0, *subDividingDimension = 0;
	    QuadDebug(Device.AllocateMemory((void **)&activeRegions, sizeof(int)*numRegions));      
	    QuadDebug(Device.AllocateMemory((void **)&subDividingDimension, sizeof(int)*numRegions));

	    if(VERBOSE)
	      {
		Println(log, "\nEntering function IntegrateFirstPhase \n");
		sprintf(msg, "\t# of input intervals\t\t: %ld\n\t#. of Thread Blocks\t\t: %ld\n\t#. of Threads per Blocks\t: %ld\n",numRegions, numBlocks, numThreads);
		Println(log, msg);
	      }



	    INTEGRATE_GPU_PHASE1<T><<<numBlocks, numThreads>>>(dRegions, dRegionsLength, numRegions, dRegionsIntegral, dRegionsError, activeRegions, subDividingDimension, epsrel, epsabs, fEvalPerRegion ,NDIM, dtest, vol);
		  
	    nregions 	+= numRegions;
	    neval 	+= numRegions*fEvalPerRegion;
		  
	    //printf("After PHASE 1 Kernel 	nregions:%i 	NumRegions:%i\n", nregions, numRegions);
	    //printf("integral:%f\n", integral);
		  
	    thrust::device_ptr<T> wrapped_ptr;
	    wrapped_ptr = thrust::device_pointer_cast(dRegionsIntegral + numRegions);
	    T rG = integral + thrust::reduce(wrapped_ptr, wrapped_ptr+numRegions);
		  
	    wrapped_ptr = thrust::device_pointer_cast(dRegionsError + numRegions);
	    T errG = error + thrust::reduce(wrapped_ptr, wrapped_ptr+numRegions);
			
	    wrapped_ptr = thrust::device_pointer_cast(dRegionsIntegral);
	    integral = integral + thrust::reduce(wrapped_ptr, wrapped_ptr+numRegions);
		  
	    wrapped_ptr = thrust::device_pointer_cast(dRegionsError);
	    error = error + thrust::reduce(wrapped_ptr, wrapped_ptr+numRegions);
		  
	    //printf("rG integral:%f integral:%f\n", rG, integral);
	    //std::cout << "Error " << errG << " " << rG<< std::endl;
	    //printf("==============================\n");
		  
	    if((errG <= MaxErr(rG, epsrel, epsabs)) && GLOBAL_ERROR) 
	      {
		if(VERBOSE)
		  {
		    sprintf(msg, "Global Error Check -\t%ld integrand evaluations so far\n%lf +- %lf ", neval, rG, errG);
		    Println(log, msg);
		  }
			
		integral = rG;
		error = errG;
		numRegions = 0;
		return;
	      }
		  
	    GenerateActiveIntervals(activeRegions, subDividingDimension);
	    QuadDebug(cudaFree(subDividingDimension));
	    QuadDebug(cudaFree(activeRegions));
	    QuadDebug(cudaFree(dRegionsError));
	    QuadDebug(cudaFree(dRegionsIntegral));      
	  }

	  template<typename integT>		
	  void IntegrateFirstPhase(T epsrel, T epsabs, T &integral, T &error, size_t &nregions, size_t &neval, integT *test, Volume *vol = NULL)
	  {
	    //Test<GPU> *dTest;
	    integT *dTest;
	    Volume *dvol;
	    if(vol)
	    {
	      cudaMalloc((void**)&dvol, sizeof(Volume));
	      cudaMemcpy(dvol, vol, sizeof(Volume), cudaMemcpyHostToDevice);
	    }
	    //x.Initialize(cInterpR, cInterpC, tau_arr, 22, 5);
	    //cudaMalloc((void**)&dTest, sizeof(Test<GPU>));
	    //cudaMemcpy(dTest, test, sizeof(Test<GPU>), cudaMemcpyHostToDevice);
	    cudaMalloc((void**)&dTest, sizeof(integT));
	    cudaMemcpy(dTest, test, sizeof(integT), cudaMemcpyHostToDevice);
	    
	    
	    //kernel<<<1,1>>>(x.interpT, x.interpT);
	    
	    cudaDeviceSynchronize();
	    for(int i  = 0; i < 100; i++){
	      FirstPhaseIteration(epsrel, epsabs, integral, error, nregions, neval, dTest, dvol);
	      if(VERBOSE){
		sprintf(msg, "Iterations %d:\t%ld integrand evaluations so far\n%lf +- %lf ", i+1 , neval, integral, error);
		Println(log, msg);
		sprintf(msg, "\n==========================================================================\n");
		Println(log, msg);
	      }
				
	      if(numRegions < 1) 
		return;
	      if(numRegions >= FIRST_PHASE_MAXREGIONS) 
		break;
	    }
	    //bring region info back to host
		  
	    hRegions = (T *)Host.AllocateMemory(&hRegions, sizeof(T) * numRegions * NDIM);
	    hRegionsLength = (T *)Host.AllocateMemory(&hRegionsLength, sizeof(T) * numRegions * NDIM);
	    QuadDebug(cudaMemcpy(hRegions, dRegions, sizeof(T) * numRegions * NDIM, cudaMemcpyDeviceToHost));
	    QuadDebug(cudaMemcpy(hRegionsLength, dRegionsLength, sizeof(T) * numRegions * NDIM, cudaMemcpyDeviceToHost));
	  }
		

	  template<typename integT>	  
	  int IntegrateSecondPhase(T epsrel, T epsabs, T &integral, T &error, size_t &nregions, size_t &neval, integT *test, T *optionalInfo = 0){
			
	    int numFailedRegions = 0;
	    int num_gpus = 0;	// number of CUDA GPUs
	    if(optionalInfo!=0){
	      optionalInfo[0] = -INFTY;
	    }
		  
	    /////////////////////////////////////////////////////////////////
	    // determine the number of CUDA capable GPUs
	    //
	    cudaGetDeviceCount(&num_gpus);
	    if(num_gpus < 1){
	      fprintf(stderr,"no CUDA capable devices were detected\n");
	      exit(1);
	    }
	    int num_cpu_procs = omp_get_num_procs();


	    /*
	      Why did you have this section?
	      for(int i = 1; i < num_gpus; i++){
	      int gpu_id;
	      QuadDebug(cudaSetDevice(i));	// "% num_gpus" allows more CPU threads than GPU devices
	      QuadDebug(cudaGetDevice(&gpu_id));
	      QuadDebug(cudaDeviceReset());
	      }	  
	    */

	    if(VERBOSE){
	      /////////////////////////////////////////////////////////////////
	      // display CPU and GPU configuration
	      sprintf(msg, "number of host CPUs:\t%d\n", omp_get_num_procs());
	      Println(log, msg);
	      sprintf(msg,"number of CUDA devices:\t%d\n", num_gpus);
	      Println(log, msg);
	      for(int i = 0; i < num_gpus; i++){
		cudaDeviceProp dprop;
		cudaGetDeviceProperties(&dprop, i);
		sprintf(msg,"   %d: %s\n", i, dprop.name);
		Println(log, msg);
	      }
	      Println(log, "---------------------------\n");
	    }

	    if(NUM_DEVICES > num_gpus)
	      NUM_DEVICES = num_gpus;
			
	    omp_set_num_threads(NUM_DEVICES);
	    cudaStream_t stream[NUM_DEVICES];
	    cudaEvent_t event[NUM_DEVICES];
			
#pragma omp parallel
	    {	
	      unsigned int cpu_thread_id 	 = omp_get_thread_num();
	      unsigned int num_cpu_threads = omp_get_num_threads();
	      if(cpu_thread_id == 0)
		printf("Phase 2 # of CPU threads:%i\n", num_cpu_threads);
	      // set and check the CUDA device for this CPU thread
	      int gpu_id = -1;
			
	      QuadDebug(cudaSetDevice(cpu_thread_id % num_gpus));	// "% num_gpus" allows more CPU threads than GPU devices
	      QuadDebug(cudaGetDevice(&gpu_id));
	      warmUpKernel<<<FIRST_PHASE_MAXREGIONS, BLOCK_SIZE>>>();


	      if(VERBOSE){
		sprintf(msg, "CPU thread %d (of %d) uses CUDA device %d\n", cpu_thread_id, num_cpu_threads, gpu_id);
		Println(log, msg);
	      }
			
	      if(cpu_thread_id < num_cpu_threads){	
		//divide regions to CPU threads
		size_t numRegionsThread = numRegions/num_cpu_threads;
				
		int startIndex = cpu_thread_id * numRegionsThread;
		int endIndex = (cpu_thread_id+1) * numRegionsThread;
				
		if(cpu_thread_id == (num_cpu_threads-1)){
		  endIndex = numRegions;
		}
				
		numRegionsThread = endIndex - startIndex;
		//QuadDebug(Device.SetHeapSize());
		CudaCheckError();
				
		Rule.loadDeviceConstantMemory(cpu_thread_id);
				
		size_t numThreads = BLOCK_SIZE;
		size_t numBlocks = numRegionsThread;
		T *dRegionsError = 0, 	*dRegionsIntegral = 0;
		T *dRegionsThread = 0, 	*dRegionsLengthThread = 0;
				
		QuadDebug(Device.AllocateMemory((void **)&dRegionsIntegral, sizeof(T)*numRegionsThread));
		QuadDebug(Device.AllocateMemory((void **)&dRegionsError, 	sizeof(T)*numRegionsThread));
		  
		int *activeRegions = 0, *subDividingDimension = 0, *dRegionsNumRegion = 0;
		QuadDebug(Device.AllocateMemory((void **)&activeRegions, sizeof(int)*numRegionsThread));      
		QuadDebug(Device.AllocateMemory((void **)&subDividingDimension, sizeof(int)*numRegionsThread));
		QuadDebug(Device.AllocateMemory((void **)&dRegionsNumRegion, sizeof(int)*numRegionsThread));      
	 

		QuadDebug(Device.AllocateMemory((void **)&dRegionsThread, sizeof(T)*numRegionsThread * NDIM));
		QuadDebug(Device.AllocateMemory((void **)&dRegionsLengthThread, sizeof(T)*numRegionsThread * NDIM));
		
		//NOTE:Copy order is important
		for(int dim = 0; dim < NDIM; ++dim){
		  QuadDebug(cudaMemcpy(dRegionsThread + dim * numRegionsThread, hRegions + dim * numRegions + startIndex, sizeof(T) * numRegionsThread, cudaMemcpyHostToDevice));
		  QuadDebug(cudaMemcpy(dRegionsLengthThread + dim * numRegionsThread, hRegionsLength + dim * numRegions + startIndex, sizeof(T) * numRegionsThread, cudaMemcpyHostToDevice));
		}

		cudaEvent_t start;
		QuadDebug(cudaStreamCreate(&stream[gpu_id]));
		QuadDebug(cudaEventCreate(&start));
		QuadDebug(cudaEventCreate(&event[gpu_id]));
		QuadDebug(cudaEventRecord(start, stream[gpu_id]));
		  
		if(VERBOSE){
		  Println(log, "\n GPU Function PHASE2");
		  sprintf(msg, "\t# of input intervals\t\t: %ld\n\t#. of Thread Blocks\t\t: %ld\n\t#. of Threads per Blocks\t: %ld\n",numRegionsThread, numBlocks, numThreads);
		  Println(log, msg);
		}
			
		//std::cout << " phase2  # of BLOCKS:" << numBlocks << " # of THREADS:" << numThreads << " " << std::endl; 
		BLOCK_INTEGRATE_GPU_PHASE2<T><<<numBlocks, numThreads, 0, stream[gpu_id]>>>(dRegionsThread, dRegionsLengthThread, numRegionsThread, dRegionsIntegral, dRegionsError, dRegionsNumRegion, activeRegions, subDividingDimension, epsrel, epsabs, fEvalPerRegion, NDIM, test);
				
		CudaCheckError();
		cudaDeviceSynchronize();
		cudaEventRecord( event[gpu_id], stream[gpu_id]);
		cudaEventSynchronize( event[gpu_id] );
		  
		float elapsed_time;
		cudaEventElapsedTime(&elapsed_time, start, event[gpu_id]);
		if(optionalInfo!=0 && elapsed_time > optionalInfo[0]){
		  optionalInfo[0] = elapsed_time;
		}

		if(VERBOSE){
		  sprintf(msg, "\nSecond Phase Kernel by thread %d (of %d) using CUDA device %d took %.1f ms ", cpu_thread_id, num_cpu_threads, gpu_id, elapsed_time);
		  Println(log, msg);
		}
			
		cudaEventDestroy(start);
		cudaEventDestroy(event[gpu_id]);

		thrust::device_ptr<T> wrapped_ptr;
		wrapped_ptr = thrust::device_pointer_cast(dRegionsIntegral);
		T integResult = thrust::reduce(wrapped_ptr, wrapped_ptr+numRegionsThread);
		  
		integral += integResult; 

		wrapped_ptr = thrust::device_pointer_cast(dRegionsError);
		error = error + thrust::reduce(wrapped_ptr, wrapped_ptr+numRegionsThread);
		  
		thrust::device_ptr<int> int_ptr = thrust::device_pointer_cast(dRegionsNumRegion);
		int regionCnt = thrust::reduce(int_ptr, int_ptr+numRegionsThread);
		nregions += regionCnt;
		//std::cout << "Num regions : " << regionCnt << std::endl;
		  
		neval += (regionCnt - numRegionsThread)*fEvalPerRegion*2+numRegionsThread*fEvalPerRegion;
	 
		int_ptr = thrust::device_pointer_cast(activeRegions);
		numFailedRegions += thrust::reduce(int_ptr, int_ptr+numRegionsThread);

		std::cout << "--" << numFailedRegions << std::endl;
		//QuadDebug(cudaThreadExit());
		  
		QuadDebug(Device.ReleaseMemory(dRegionsError));
		QuadDebug(Device.ReleaseMemory(dRegionsIntegral));
		QuadDebug(Device.ReleaseMemory(dRegionsThread));
		QuadDebug(Device.ReleaseMemory(dRegionsLengthThread));
		QuadDebug(Device.ReleaseMemory(activeRegions));
		QuadDebug(Device.ReleaseMemory(subDividingDimension));
		QuadDebug(Device.ReleaseMemory(dRegionsNumRegion));
		QuadDebug(cudaDeviceSynchronize());
	      }	
	    }
		
	    //sprintf(msg, "Execution time : %.2lf", optionalInfo[0]);
	    //Print(msg);
	    return numFailedRegions;

	  }

	};
   
}
