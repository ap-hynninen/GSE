#ifndef CPUGRID_H
#define CPUGRID_H
//
// CPU grid class
// (c) Antti-Pekka Hynninen Feb 2015
// aphynninen@hotmail.com
//
#include <cmath>
#include "cuda_utils.h"
#include "Grid.h"

// Forward declaration of CudaGrid
template<typename T> class CudaGrid;

template <typename T>
class CpuGrid : public Grid<T> {
private:

public:
  CpuGrid(const int sizeX, const int sizeY, const int sizeZ) :
    Grid<T>(sizeX, sizeY, sizeZ) { this->data = new T[this->size]; }

  ~CpuGrid() { delete [] this->data; }

  void clear() {
    for (int i=0;i < this->size;i++) this->data[i] = (T)0;
  }

  void add(const Grid<T>& grid) {
    const T* gridData = grid.getDataPointer();
    for (int i=0;i < this->size;i++) this->data[i] += gridData[i];
  }

  void sub(const Grid<T>& grid) {
    const T* gridData = grid.getDataPointer();
    for (int i=0;i < this->size;i++) this->data[i] -= gridData[i];
  }

  void copy(const CpuGrid<T>& grid) {
    const T* gridData = grid.getDataPointer();
    for (int i=0;i < this->size;i++) this->data[i] = gridData[i];
  }
  
  void copy(const CudaGrid<T>& d_grid) {
    const T* d_gridData = d_grid.getDataPointer();
    copy_DtoH_sync<T>(d_gridData, this->data, this->size);
  }

  
  void scale(const T val) {
    for (int i=0;i < this->size;i++) this->data[i] *= val;
  }

  T sum() {
    T res = (T)0.0;
    for (int i=0;i < this->size;i++) res += this->data[i];
    return res;
  }
  
  T maxAbsValue() {
    T maxval = (T)0;
    for (int i=0;i < this->size;i++) {
      T val = fabs(this->data[i]);
      maxval = (maxval < val) ? val : maxval;
    }
    return maxval;
  }

  // Get value to CPU
  T getDataValue(const int x, const int y, const int z) {return this->data[this->getPos(x,y,z)];}
  T getDataValue(const int x, const int y, const int z) const {return this->data[this->getPos(x,y,z)];}

  void save(const char* filename) const {
    std::ofstream file(filename);
    if (file.is_open()) {
      for (int z=0;z < this->sizeZ;z++) {
	for (int y=0;y < this->sizeY;y++) {
	  for (int x=0;x < this->sizeX;x++) {
	    file << getDataValue(x, y, z) << std::endl;
	  }
	}
      }
    } else {
      std::cerr << "CpuGrid<T>::save, Error opening file " << filename << std::endl;
      exit(1);
    }
  }

  void save(const char* filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
      for (int z=0;z < this->sizeZ;z++) {
	for (int y=0;y < this->sizeY;y++) {
	  for (int x=0;x < this->sizeX;x++) {
	    file << getDataValue(x, y, z) << std::endl;
	  }
	}
      }
    } else {
      std::cerr << "CpuGrid<T>::save, Error opening file " << filename << std::endl;
      exit(1);
    }
  }

};

#endif // CPUGRID_H
