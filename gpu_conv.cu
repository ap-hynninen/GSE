#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include "TypesGSE.h"
#include "cuda/cuda_utils.h"
#include "CudaConvolution.h"
#include "CpuConvolution.h"
#include "CpuGrid.h"
#include "CudaGrid.h"

template<typename T>
void test(const int n, const int kernelRadius);

int main() {
    
  std::vector<int> devices;
  start_gpu(1, 0, devices);

  test<double>(128, 5);
  
  stop_gpu();
  return 0;
}

template<typename T>
void test(const int n, const int kernelRadius) {
  T* kernelX = new T[2*kernelRadius+1];
  T* kernelY = new T[2*kernelRadius+1];
  T* kernelZ = new T[2*kernelRadius+1];

  const double sigma = 2.12;
  const double inv_2sigmasq = 1.0/(2.0*sigma*sigma);
  const double pref = pow(2*pi_dbl*sigma*sigma,-3.0/2.0);
  
  for (int ii=-kernelRadius;ii <= kernelRadius;ii++) {
    kernelX[ii+kernelRadius] = exp(-(ii*ii)*inv_2sigmasq);
  }
  for (int ii=-kernelRadius;ii <= kernelRadius;ii++) {
    kernelY[ii+kernelRadius] = exp(-(ii*ii)*inv_2sigmasq);
  }
  for (int ii=-kernelRadius;ii <= kernelRadius;ii++) {
    kernelZ[ii+kernelRadius] = pref*exp(-(ii*ii)*inv_2sigmasq);
  }
  
  CpuGrid<T> cpuRho(n, n, n);
  CpuGrid<T> cpuRhoXYZ(n, n, n);
  CudaGrid<T> cudaRho(n, n, n);
  CudaGrid<T> cudaRhoXYZ(n, n, n);

  srand(12301);
  T* cpuRhoData = cpuRho.getDataPointer();
  for (int i=0;i < n*n*n;i++) {
    cpuRhoData[i] = (T)((double)rand()/((double)RAND_MAX + 1.0));
  }
  cudaRho.copy(cpuRho);

  CudaConvolution<T> cudaConv(n, n, n);
  cudaConv.setupKernel(5, 5, 5, kernelX, kernelY, kernelZ);
  cudaConv.run(cudaRho, cudaRhoXYZ);

  CpuConvolution<T> cpuConv(n, n, n);
  cpuConv.setupKernel(5, 5, 5, kernelX, kernelY, kernelZ);
  cpuConv.run(cpuRho, cpuRhoXYZ);

  CpuGrid<T> cpuDiffXYZ(n, n, n);  
  cpuDiffXYZ.copy(cudaRhoXYZ);
  cpuDiffXYZ.sub(cpuRhoXYZ);
  T maxDiff = cpuDiffXYZ.maxAbsValue();

  printf("maxDiff = %lf\n",maxDiff);
  
  delete [] kernelX;
  delete [] kernelY;
  delete [] kernelZ;
}
