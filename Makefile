NVCC=/usr/local/cuda-10.0/bin/nvcc -x cu
CCFLAGS=-arch=sm_70 \
        -std=c++14 \
        -I /usr/include/mpich-ppc64le \
        -I /opt/gsl-2.5/include \
        -Xcompiler -fopenmp

LDFLAGS=-L/usr/local/gcc-6.4.0/lib -L/usr/lib64/mpich/lib -lmpich -lgomp

all: cudaCuhre toy

toy: toy.cu
	$(NVCC) $(CCFLAGS) -o $@ $< $(LDFLAGS)

cudaCuhre: cudaCuhre.cu lc_lt_t.h quad/GPUquad/GPUquad.h
	$(NVCC) $(CCFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f cudaCuhre toy *.o *.log

run:
	./cudaCuhre

