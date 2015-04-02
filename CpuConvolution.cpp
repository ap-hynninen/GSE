
#include <cassert>
#include "CpuConvolution.h"

//
// Class creator
//
template <typename T>
CpuConvolution<T>::CpuConvolution(const int sizeX, const int sizeY, const int sizeZ) :
  dataX(sizeX, sizeY, sizeZ), dataXY(sizeX, sizeY, sizeZ) {

  kernelX = NULL;
  kernelXlen = 0;
  kernelY = NULL;
  kernelYlen = 0;
  kernelZ = NULL;
  kernelZlen = 0;
}

//
// Class destructor
//
template <typename T>
CpuConvolution<T>::~CpuConvolution() {
  if (kernelX != NULL) delete [] kernelX;
  if (kernelY != NULL) delete [] kernelY;
  if (kernelZ != NULL) delete [] kernelZ;
}

//
// Setup kernel
//
template <typename T>
void CpuConvolution<T>::setupKernel(const int kernelRadiusX, const int kernelRadiusY, const int kernelRadiusZ,
				     const T* h_kernelX, const T* h_kernelY, const T* h_kernelZ) {
  this->kernelRadiusX = kernelRadiusX;
  this->kernelRadiusY = kernelRadiusY;
  this->kernelRadiusZ = kernelRadiusZ;

  if (kernelX != NULL && kernelXlen < 2*kernelRadiusX+1) {
    delete [] kernelX;
    kernelX = NULL;
  }
  if (kernelX == NULL) {
    kernelXlen = 2*kernelRadiusX+1;
    kernelX = new T[kernelXlen];
  }

  if (kernelY != NULL && kernelYlen < 2*kernelRadiusY+1) {
    delete [] kernelY;
    kernelY = NULL;
  }
  if (kernelY == NULL) {
    kernelYlen = 2*kernelRadiusY+1;
    kernelY = new T[kernelYlen];
  }

  if (kernelZ != NULL && kernelZlen < 2*kernelRadiusZ+1) {
    delete [] kernelZ;
    kernelZ = NULL;
  }
  if (kernelZ == NULL) {
    kernelZlen = 2*kernelRadiusZ+1;
    kernelZ = new T[kernelZlen];
  }

  for (int i=0;i < 2*kernelRadiusX+1;i++) kernelX[i] = h_kernelX[i];
  for (int i=0;i < 2*kernelRadiusY+1;i++) kernelY[i] = h_kernelY[i];
  for (int i=0;i < 2*kernelRadiusZ+1;i++) kernelZ[i] = h_kernelZ[i];
}

//
// Run convolution
//
template <typename T>
void CpuConvolution<T>::run(const CpuGrid<T>& data, CpuGrid<T>& dataXYZ) {  
  const int sizeX = dataX.getSizeX();
  const int sizeY = dataX.getSizeY();
  const int sizeZ = dataX.getSizeZ();
  assert(sizeX == data.getSizeX());
  assert(sizeY == data.getSizeY());
  assert(sizeZ == data.getSizeZ());
  assert(sizeX == dataXYZ.getSizeX());
  assert(sizeY == dataXYZ.getSizeY());
  assert(sizeZ == dataXYZ.getSizeZ());

  dataX.clear();
  dataXY.clear();
  dataXYZ.clear();

  const T* datap = data.getDataPointer();
  T* dataXp = dataX.getDataPointer();

  // Convolute X
  for (int pz=0;pz < sizeZ;pz++) {
    for (int py=0;py < sizeY;py++) {
      int pyz = (py + pz*sizeY)*sizeX;
      for (int ix=0;ix < sizeX;ix++) {
	T sum = 0;
	for (int ii=-kernelRadiusX;ii <= kernelRadiusX;ii++) {
	  int px = (ii+ix+sizeX) % sizeX;
	  int p = px + pyz;
	  sum += kernelX[ii+kernelRadiusX]*datap[p];
	}
	dataXp[ix + pyz] = sum;
      }
    }
  }

  dataXYZ.copy(dataX);
}

//
// Explicit instances of this class
//
template class CpuConvolution<float>;
template class CpuConvolution<double>;
