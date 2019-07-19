# Add source files here
EXECUTABLE      := output 
# Cuda source files 
CUFILES         := cudaCuhre.cu integrands/interp_2d.cc

NVCC            := /usr/local/cuda-10.0/bin/nvcc
CCFLAGS         := -lgomp -Xcompiler -fopenmp -arch=sm_70 -std=c++14
INCLUDE  	:= -I /usr/include/mpich-ppc64le/ -I include  
LIBS     	:= -L/usr/local/gcc-6.4.0/lib -L/usr/lib64/mpich/lib -lmpich -L lib/ -lgsl -lgslcblas 

export LD_LIBRARY_PATH=/usr/local/gcc-6.4.0/lib

${EXECUTABLE}: ${CUFILES}
	${NVCC} -x cu -o ${EXECUTABLE} ${INCLUDE} ${LIBS} ${CCFLAGS} ${CUFILES}
clean:
	-$(RM)	${EXECUTABLE}\
	*.o *.*log
run:
	./output

