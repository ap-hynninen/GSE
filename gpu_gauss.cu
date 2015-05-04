#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include "TypesGSE.h"
#include "cuda/cuda_utils.h"
#include "CpuGrid.h"
#include "CudaGrid.h"
#include "CpuGaussCharge.h"
#include "CudaGaussCharge.h"

template<typename T>
void test();

int main() {
  std::vector<int> devices;
  start_gpu(1, 0, devices);

  test<float>();
  
  stop_gpu();
  return 0;
}

//
// Loads array data from file
//
template <typename T>
void loadArray(const int width, const int numLine, const char *filename, T *array) {
  std::ifstream file(filename);
  if (file.is_open()) {

    for (int i=0;i < numLine;i++) {
      for (int k=0;k < width;k++) {
	if (!(file >> array[i*width+k])) {
	  std::cerr<<"Error reading file "<<filename<<std::endl;
	  exit(1);
	}
      }
    }

  } else {
    std::cerr<<"Error opening file "<<filename<<std::endl;
    exit(1);
  }
}

//
//
//
template<typename T>
void test() {
  const double sigma = 2.12132034355964/sqrt(2.0);
  const double lambdaSigma = 3.0;
  const int n = 64;
  const int numCoord = 23558;
  const double boxx = 62.23;
  const double boxy = 62.23;
  const double boxz = 62.23;
  
  xyzq_t<T>* h_xyzq = new xyzq_t<T>[numCoord];
  loadArray<T>(4, numCoord, "data/xyzq.txt", (T *)h_xyzq);

  xyzq_t<double>* h_xyzqDBL = new xyzq_t<double>[numCoord];
  loadArray<double>(4, numCoord, "data/xyzq.txt", (double *)h_xyzqDBL);

  xyzq_t<T>* d_xyzq = NULL;
  allocate< xyzq_t<T> >(&d_xyzq, numCoord);
  copy_HtoD_sync< xyzq_t<T> >(h_xyzq, d_xyzq, numCoord);
  
  CpuGrid<T> cpuRho(n, n, n);
  CpuGrid<double> cpuRhoDBL(n, n, n);
  CudaGrid<T> cudaRho(n, n, n);

  CpuGaussCharge<T,T> cpuGaussCharge(n, n, n);
  CpuGaussCharge<double,double> cpuGaussChargeDBL(n, n, n);
  CudaGaussCharge<T,T> cudaGaussCharge(n, n, n);

  int nn = numCoord;
  cpuGaussCharge.spreadChargeToGrid(sigma, lambdaSigma*sqrt(2.0)*sigma, nn, h_xyzq,
				    boxx, boxy, boxz, cpuRho);

  cpuGaussChargeDBL.spreadChargeToGrid(sigma, lambdaSigma*sqrt(2.0)*sigma, nn, h_xyzqDBL,
				       boxx, boxy, boxz, cpuRhoDBL);

  cudaGaussCharge.spreadChargeToGrid(sigma, lambdaSigma*sqrt(2.0)*sigma, nn, d_xyzq,
				     boxx, boxy, boxz, cudaRho);

  CpuGrid<T> cpuRhoTmp(n, n, n);

  cpuRhoTmp.copy(cpuRhoDBL);
  cpuRhoTmp.sub(cpuRho);
  T maxDiff = cpuRhoTmp.maxAbsValue();
  printf("max|cpuRho - cpuRhoDBL| = %e\n",maxDiff);
  
  cpuRhoTmp.copy(cudaRho);
  cpuRhoTmp.sub(cpuRho);
  maxDiff = cpuRhoTmp.maxAbsValue();
  printf("max|cpuRho - cudaRho| = %e\n",maxDiff);

  cpuRho.save("cpuRho.txt");
  cpuRhoTmp.copy(cudaRho);
  cpuRhoTmp.save("cudaRho.txt");
  
  delete [] h_xyzq;
  delete [] h_xyzqDBL;
  deallocate< xyzq_t<T> >(&d_xyzq);
}
