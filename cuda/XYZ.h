#ifndef XYZ_H
#define XYZ_H

#include <iostream>
#include <fstream>

//
// XYZ array base class
//
// (c) Antti-Pekka Hynninen, 2014, aphynninen@hotmail.com
//

template <typename T>
class XYZ {

protected:

  int _size;     // Size of the XYZ arrays
  int _capacity; // Capacity of the XYZ arrays
  T* _x;         // X-array
  T* _y;         // Y-array
  T* _z;         // Z-array

 public:
  
  // Default constructor
  XYZ() {
    _size = 0;
    _capacity = 0;
    _x = NULL;
    _y = NULL;
    _z = NULL;
  }

  // "Create from existing pointer" -constructor
  XYZ(int size, int capacity, T *x, T *y, T *z) {
    this->_size = _size;
    this->_capacity = _capacity;
    this->_x = x;
    this->_y = y;
    this->_z = z;
  }

  // Returns true if the XYZ strided arrays match in data content sizes
  template <typename P>
  bool match(const XYZ<P> &xyz) {
    return ((sizeof(T) == sizeof(P)) && (this->_size == xyz._size));
  }

  // Swaps XYZ contents
  void swap(XYZ<T> &xyz) {
    assert(this->match(xyz));
    
    // Swap pointers
    std::swap(this->_x, xyz._x);
    std::swap(this->_y, xyz._y);
    std::swap(this->_z, xyz._z);

    // Swap size & capacity
    std::swap(this->_size, xyz._size);
    std::swap(this->_capacity, xyz._capacity);
  }

  // capacity = current capacity
  // size     = desired size
  virtual void realloc_array(T** array, int* capacity, int size, float fac) = 0;
  virtual void resize_array(T** array, int* capacity, int size, int new_size, float fac) = 0;

  // Resize array, preserves content
  void resize(int size, float fac=1.0f) {
    int x_capacity = this->_capacity;
    int y_capacity = this->_capacity;
    int z_capacity = this->_capacity;
    resize_array(&this->_x, &x_capacity, this->_size, size, fac);
    resize_array(&this->_y, &y_capacity, this->_size, size, fac);
    resize_array(&this->_z, &z_capacity, this->_size, size, fac);
    assert(x_capacity == y_capacity);
    assert(x_capacity == z_capacity);
    this->_size = size;
    this->_capacity = x_capacity;
  }

  // Re-allocates array, does not preserve content
  void realloc(int size, float fac=1.0f) {
    int x_capacity = this->_capacity;
    int y_capacity = this->_capacity;
    int z_capacity = this->_capacity;
    realloc_array(&this->_x, &x_capacity, size, fac);
    realloc_array(&this->_y, &y_capacity, size, fac);
    realloc_array(&this->_z, &z_capacity, size, fac);
    assert(x_capacity == y_capacity);
    assert(x_capacity == z_capacity);
    this->_size = size;
    this->_capacity = x_capacity;
  }

  virtual void get_host_xyz(T*& hx, T*& hy, T*& hz)=0;
  virtual void release_host_xyz(T*& hx, T*& hy, T*& hz)=0;
  
  // Save to file
  void save(const char* filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
      T *hx, *hy, *hz;
      get_host_xyz(hx, hy, hz);
      for (int i=0;i < _size;i++) {
	file << hx[i] << " " << hy[i] << " " << hz[i] << std::endl;
      }
      release_host_xyz(hx, hy, hz);
    } else {
      std::cerr<<"Error opening file "<<filename<<std::endl;
      exit(1);
    }
  }

  int size() {return _size;}
  int size() const {return _size;}

  T* x() {return _x;}
  T* y() {return _y;}
  T* z() {return _z;}
  const T* x() const {return _x;}
  const T* y() const {return _y;}
  const T* z() const {return _z;}

};

#endif // XYZ_H
