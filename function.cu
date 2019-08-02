/*
#define FUN 5
#define DIM 8

#ifndef FUN
        #define FUN 5
        #define DIM 4
#endif
*/

#define PI 3.14159265358979323844
#define MIN(a, b) (((a) < (b)) ? a : b)

//#define TYPE double
__device__ float Interp(double y,
                        double x,
                        cudaTextureObject_t interpT,
                        cudaTextureObject_t interpR,
                        cudaTextureObject_t interpC,
                        int offset);
__host__ float Interp(double y,
                      double x,
                      float* interpT,
                      float* interpR,
                      float* interpC);

typedef bool GPU;
typedef int CPU;

template <typename T>
__device__ T
r8vec_sum(int n, const T a[])
{
  T sum;
  sum = 0.0;
  for (int i = 0; i < n; i++) {
    sum = sum + a[i];
  }
  return sum;
}

template <typename T>
__device__ T
r8vec_dot(int n, T a1[], const T a2[])
{
  int i;
  T value;

  value = 0.0;
  for (i = 0; i < n; i++) {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//==================================
// added by Ioannis
/* __device__ double gaussian(double x, double mu, double sigma)
{
    double z = (x - mu) / sigma;
    return exp(-z * z / 2.) * 0.3989422804014327 / sigma;
}


__device__ double erfc_scaled(double a, double b, double root_two_sigma){
        return  erfc((a - b) / root_two_sigma);
        }*/

//==================================

template <typename T>
__device__ T
IntegrandFunc(const T xx[], int NDIM)
{
  T f = 0;

#if FUN == 1
  T t1 = 0;
  int N = ;
  for (N = 1; N <= NDIM; ++N) {
    t1 += pow(xx[N - 1], 2);
  }
  f = 1.0 / (0.1 + pow(cos(t1), 2));

#elif FUN == 2
  // FTEST3
  T t1 = 1.0;
  int N = 1;
  for (N = 1; N <= NDIM; ++N) {
    t1 = t1 * cos(pow(2.0, 2.0 * N) * xx[N - 1]);
  }
  f = cos(t1);

#elif FUN == 3
  // FTEST6
  T t1 = 1.0;
  int N = 1;
  for (N = 1; N <= NDIM; ++N) {
    t1 = t1 * N * asin(pow(xx[N - 1], N));
  }
  f = sin(t1);

#elif FUN == 4
  // FTEST7
  T t1 = 1.0;
  int N = 1;
  for (N = 1; N <= NDIM; ++N) {
    t1 = t1 * asin(xx[N - 1]);
  }
  f = sin(t1);

#elif FUN == 5
  T sum = 0;
  int N = 1;
  for (N = 1; N <= NDIM; ++N) {
    sum = sum - cos(10.0 * xx[N - 1]) / 0.054402111088937;
  }
  f = sum / 2.0;

#elif FUN == 6
  T sum = 0;
  int N = 1;
  T alpha = 110.0 / (NDIM * NDIM * sqrt(NDIM * 1.0));
  for (N = 1; N <= NDIM; ++N) {
    sum += alpha * xx[N - 1];
  }

  sum = 2.0 * PI * alpha + sum;
  f = cos(sum);
#elif FUN == 7
  int N = 1;
  T total = 1.0;
  T alpha = 600.0 / (NDIM * NDIM * NDIM);
  T beta = alpha;
  for (N = 1; N <= NDIM; ++N) {
    total = total * (1.0 / pow(alpha, 2) + pow(xx[N - 1] - beta, 2));
  }
  f = 1.0 / total;

#elif FUN == 10
  T total = 0.0;
  T alpha = 0.5; // 150.0/(NDIM * NDIM * NDIM);
  T beta = alpha;
  int N = 1;

  /*  T tau 		=  	Interp(xx[0], xx[1], interpT, interpR, interpC,
    0); T mu 		= 	Interp(xx[3], xx[2], interpT, interpR, interpC,
    1); T sigma 		=	Interp(xx[7], xx[5], interpT, interpR,
    interpC, 2); T fmsk 		= 	Interp(xx[3], xx[6], interpT,
    interpR, interpC, 3); T fprj 		=	Interp(xx[0], xx[5],
    interpT, interpR, interpC, 4);
   */

  for (N = 1; N <= NDIM; ++N) {
    total = total + alpha * fabs(xx[N - 1] - beta);
  }
  f = exp(-total);

#elif FUN == 11

  /*T highs[3] = {2, .3, 1};
        T lows[3] =  {1, .1, 0};
        
        
        
        //adjust x values according to range
         T lt = (highs[0]-lows[0])*xx[0] + lows[0];
         T zt = (highs[1]-lows[1])*xx[1] + lows[1];
         T lc = (highs[2]-lows[2])*xx[2] + lows[2];
         
         //test interpolation values
         //lt = 2;
        // zt = .12;
        
         T tau 		=  	Interp(lt, zt, interpT, interpR, interpC, 0);
         T mu 		= 	Interp(lt, zt, interpT, interpR, interpC, 1);
     T sigma 	=	Interp(lt, zt, interpT, interpR, interpC, 2);
     T fmsk 	= 	Interp(lt, zt, interpT, interpR, interpC, 3);
     T fprj 	=	Interp(lt, zt, interpT, interpR, interpC, 4);
         
         if(threadIdx.x == 0 && blockIdx.x==0)
         {
                 printf("tau:%f mu:%f sigma:%f fmsk:%f fprj:%f\n", tau, mu,
     sigma, fmsk, fprj);
         }
                
         T exptau = exp(tau * (2.0 * mu + tau * sigma * sigma - 2.0 * lc)
     / 2.0); T root_two_sigma = sqrt(2.0) * sigma; T mu_tau_sig_sqr = mu + tau *
     sigma * sigma;
         
     f = (1.0 - fmsk) * (1.0 - fprj) * gaussian(lc, mu, sigma) +
              0.5 * ((1.0 - fmsk) * fprj * tau + fmsk * fprj / lt) * exptau *
     erfc_scaled(mu_tau_sig_sqr, lc, root_two_sigma) + 0.5 * fmsk / lt *
     (erfc_scaled(lc, mu, root_two_sigma) - erfc_scaled(lc + lt, mu,
     root_two_sigma)) - 0.5 * fmsk * fprj / lt *(exp(-tau * lt) * exptau *
     erfc_scaled(mu_tau_sig_sqr, lc + lt, root_two_sigma));
        
        for(int i=0; i<DIM; i++)
                f *= (highs[i]-lows[i]);
*/
#elif FUN == 12
  float a = 1.0;
  float b = 2.0;
  if (threadIdx.x == 0 && blockIdx.x == 0)
    printf("x: %f	%f 		%f\n", xx[0], xx[1], xx[2]);
  float x1 = (b - a) * xx[0] + a;
  float x2 = (b - a) * xx[1] + a;
  float x3 = (3 - 2) * xx[2] + 2;
  if (threadIdx.x == 0 && blockIdx.x == 0)
    printf("x: %f	%f 		%f\n", x1, x2, x3);
  f = (x1 * x2 * x3) * (b - a) * (3 - 2);
#endif

  return f;
}
