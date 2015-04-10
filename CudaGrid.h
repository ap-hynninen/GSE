#ifndef CUDAGRID_H
#define CUDAGRID_H
//
// CUDA grid class
// (c) Antti-Pekka Hynninen Feb 2015
// aphynninen@hotmail.com
//
#include "cuda/cuda_utils.h"
#include "Grid.h"

// Forward declaration of CpuGrid
template<typename T> class CpuGrid;

template <typename T>
class CudaGrid : public Grid<T> {
private:
 
public:
  CudaGrid(const int sizeX, const int sizeY, const int sizeZ) :
    Grid<T>(sizeX, sizeY, sizeZ) {
    allocate<T>(&this->data, this->size);
  }

  ~CudaGrid() {
    deallocate<T>(&this->data);
  }

  void clear() {
    clear_gpu_array<T>(this->data, this->size);
  }
  
  void copy(const CpuGrid<T>& h_grid) {
    const T* h_gridData = h_grid.getDataPointer();
    copy_HtoD_sync<T>(h_gridData, this->data, this->size);
  }

  void copy(const CudaGrid<T>& d_grid) {
    const T* d_gridData = d_grid.getDataPointer();
    copy_DtoD_sync<T>(d_gridData, this->data, this->size);
  }

};

#endif // CUDAGRID_H
