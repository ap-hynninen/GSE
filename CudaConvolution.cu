
#include <cassert>
#include <cuda.h>
#include "CudaConvolution.h"

//
// Convolve in x-direction
// data  = input data
// dataX = output data convolved in X direction
//
template <typename T>
__global__ void convXkernel(const int sizeX, const int sizeY, const int sizeZ,
			    const int numStep,
			    const int kernelRadius,
			    const T* __restrict__ kernel,
			    const T* __restrict__ data,
			    T* __restrict__ dataX) {
  // Shared memory requirement:
  // sizeof(T)*( (numStep*blockDim.x + 2*kernelRadius)*blockDim.y + 2*kernelRadius+1)
  extern __shared__ char sh_bufX[];
  T* sh_data = (T*)&sh_bufX[0];
  T* sh_kernel = (T*)&sh_bufX[(numStep*blockDim.x + 2*kernelRadius)*blockDim.y*sizeof(T)];
  
  const int baseX = blockIdx.x*numStep*blockDim.x;
  const int baseY = blockIdx.y*blockDim.y + threadIdx.y;
  const int baseZ = blockIdx.z;
  const int baseYZ = (baseY + baseZ*sizeY)*sizeX;
  const int base = baseX + baseYZ;
  const int sh_base = threadIdx.y*(numStep*blockDim.x + 2*kernelRadius);

  // Load kernel into shared memory
  for (int j=threadIdx.x+threadIdx.y*blockDim.x;j < 2*kernelRadius+1;j+=blockDim.x*blockDim.y) {
    sh_kernel[j] = kernel[j];
  }
  
  // Load main data into shared memory
  for (int i=0;i < numStep;i++) {
    sh_data[sh_base + kernelRadius + i*blockDim.x + threadIdx.x] = data[base + (i*blockDim.x + threadIdx.x)%sizeX];
  }

  // Load left kernel radius into shared memory
  for (int i=threadIdx.x;i < kernelRadius;i+=blockDim.x) {
    sh_data[sh_base + i] = data[base + ( -kernelRadius + i + sizeX)%sizeX];
  }

  // Load right kernel radius into shared memory
  for (int i=threadIdx.x;i < kernelRadius;i+=blockDim.x) {
    sh_data[sh_base + kernelRadius + numStep*blockDim.x + i] = data[base + ( numStep*blockDim.x + i)%sizeX];
  }

  // Sync
  __syncthreads();

  // Convolve and store result
  for (int i=0;i < numStep;i++) {
    T sum = 0;
    for (int j=-kernelRadius;j <= kernelRadius;j++) {
      sum += sh_kernel[j+kernelRadius]*sh_data[sh_base + threadIdx.x + i*blockDim.x + j+kernelRadius];
    }
    dataX[base + threadIdx.x + i*blockDim.x] = sum;
  }  
  
}

//#########################################################################
//#########################################################################
//#########################################################################

//
// Class creator
//
template <typename T>
CudaConvolution<T>::CudaConvolution(const int sizeX, const int sizeY, const int sizeZ) :
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
CudaConvolution<T>::~CudaConvolution() {
  if (kernelX != NULL) deallocate<T>(&kernelX);
  if (kernelY != NULL) deallocate<T>(&kernelY);
  if (kernelZ != NULL) deallocate<T>(&kernelZ);
}

//
// Setup kernel
//
template <typename T>
void CudaConvolution<T>::setupKernel(const int kernelRadiusX, const int kernelRadiusY, const int kernelRadiusZ,
				     const T* h_kernelX, const T* h_kernelY, const T* h_kernelZ) {
  this->kernelRadiusX = kernelRadiusX;
  this->kernelRadiusY = kernelRadiusY;
  this->kernelRadiusZ = kernelRadiusZ;
  reallocate<T>(&kernelX, &kernelXlen, 2*kernelRadiusX+1);
  reallocate<T>(&kernelY, &kernelYlen, 2*kernelRadiusY+1);
  reallocate<T>(&kernelZ, &kernelZlen, 2*kernelRadiusZ+1);
  copy_HtoD_sync<T>(h_kernelX, kernelX, 2*kernelRadiusX+1);
  copy_HtoD_sync<T>(h_kernelY, kernelY, 2*kernelRadiusY+1);
  copy_HtoD_sync<T>(h_kernelZ, kernelZ, 2*kernelRadiusZ+1);
}

//
// Run convolution
//
template <typename T>
void CudaConvolution<T>::run(const CudaGrid<T>& data, CudaGrid<T>& dataXYZ) {  
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
  
  int blockDimX = 16;
  int blockDimY = 4;
  int numStep = 8;

  printf("sizeX = %d\n",sizeX);
  
  dim3 nblock(sizeX / (numStep * blockDimX), sizeY / blockDimY, sizeZ);
  dim3 nthread(blockDimX, blockDimY, 1);
  int shmem_size = sizeof(T)*( (numStep*blockDimX + 2*kernelRadiusX)*blockDimY + 2*kernelRadiusX+1);

  printf("nblock = %d %d %d\n",nblock.x,nblock.y,nblock.z);
  printf("nthread = %d %d %d\n",nthread.x,nthread.y,nthread.z);
  printf("shmem_size = %d\n",shmem_size);

  convXkernel<T> <<< nblock, nthread, shmem_size >>>
    (sizeX, sizeY, sizeZ, numStep, kernelRadiusX, kernelX,
     data.getDataPointer(), dataX.getDataPointer());
  cudaCheck(cudaGetLastError());

  cudaCheck(cudaDeviceSynchronize());
  dataXYZ.copy(dataX);
}

//
// Explicit instances of this class
//
template class CudaConvolution<float>;
template class CudaConvolution<double>;
