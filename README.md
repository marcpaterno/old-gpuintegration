# gpuintegration
Numerical integration code for GPUs.

gcc prior version 7 must be placed in the path, compatible with c++14
gcc 6.4.0 is confirmed to work

export PATH=/usr/local/gcc-6.4.0/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/gcc-6.4.0/lib:$LD_LIBRARY_PATH

in the absence of modules on imbpower9, the following library paths must be set

1. MPI. export LD_LIBRARY_PATH=/usr/lib64/mpich/lib:$LD_LIBRARY_PATH



2. gls  export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH 

3. export LD_LIBRARY_PATH=/usr/local/gcc-6.4.0/lib64:$LD_LIBRARY_PATH

The necessary headers for cuba and gsl must also be added to an "include" directory