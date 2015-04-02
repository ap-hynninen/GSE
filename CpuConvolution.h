#ifndef CPUCONVOLUTION_H
#define CPUCONVOLUTION_H
//
// CPU 3d convolution for separable kernels
// (c) Antti-Pekka Hynninen
//
#include "CpuGrid.h"

template <typename T>
class CpuConvolution {
private:

  // Temporary arrays
  CpuGrid<T> dataX;
  CpuGrid<T> dataXY;

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
  CpuConvolution(const int sizeX, const int sizeY, const int sizeZ);
  ~CpuConvolution();
  void setupKernel(const int kernelRadiusX, const int kernelRadiusY, const int kernelRadiusZ,
		   const T* h_kernelX, const T* h_kernelY, const T* h_kernelZ);
  void run(const CpuGrid<T>& data, CpuGrid<T>& dataXYZ);
};

#endif //CPUCONVOLUTION_H
