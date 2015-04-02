#ifndef CUDACONVOLUTION_H
#define CUDACONVOLUTION_H
//
// CUDA 3d convolution for separable kernels
// (c) Antti-Pekka Hynninen
//
#include "CudaGrid.h"

template <typename T>
class CudaConvolution {
private:

  // Temporary arrays
  CudaGrid<T> dataX;
  CudaGrid<T> dataXY;

  // Kernel arrays
  int kernelXlen;
  T* kernelX;
  
  int kernelYlen;
  T* kernelY;

  int kernelZlen;
  T* kernelZ;

  // Kernel radiuses
  int kernelRadiusX;
  int kernelRadiusY;
  int kernelRadiusZ;
  
public:
  CudaConvolution(const int sizeX, const int sizeY, const int sizeZ);
  ~CudaConvolution();
  void setupKernel(const int kernelRadiusX, const int kernelRadiusY, const int kernelRadiusZ,
		   const T* h_kernelX, const T* h_kernelY, const T* h_kernelZ);
  void run(const CudaGrid<T>& data, CudaGrid<T>& dataXYZ);
};

#endif //CUDACONVOLUTION_H
