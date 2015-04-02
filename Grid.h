#ifndef GRID_H
#define GRID_H
//
// Grid base class
// (c) Antti-Pekka Hynninen Feb 2015
// aphynninen@hotmail.com
//
#include <iostream>
#include <fstream>

template <typename T>
class Grid {
  ///private:
protected:
  // Size of the grid
  const int sizeX;
  const int sizeY;
  const int sizeZ;
  const int size;

  // Grid data
  T* data;

public:
  Grid(const int sizeX, const int sizeY, const int sizeZ) :
    sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), size(sizeX*sizeY*sizeZ) {}

  T* getDataPointer() {return data;}
  const T* getDataPointer() const {return data;}

  // Get safe position
  inline int getPos(const int x, const int y, const int z) {
    const int xt = (x + sizeX) % sizeX;
    const int yt = (y + sizeY) % sizeY;
    const int zt = (z + sizeZ) % sizeZ;
    return xt + (yt + zt*sizeY)*sizeX;
  }

  inline int getPos (const int x, const int y, const int z) const {
    const int xt = (x + sizeX) % sizeX;
    const int yt = (y + sizeY) % sizeY;
    const int zt = (z + sizeZ) % sizeZ;
    return xt + (yt + zt*sizeY)*sizeX;
  }
    
  int getSizeX() {return sizeX;}
  int getSizeY() {return sizeY;}
  int getSizeZ() {return sizeZ;}
  int getSize() {return size;}  
  int getSizeX() const {return sizeX;}
  int getSizeY() const {return sizeY;}
  int getSizeZ() const {return sizeZ;}
  int getSize() const {return size;}

};

#endif // GRID_H
