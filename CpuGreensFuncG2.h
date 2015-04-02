#ifndef CPUGREENSFUNCG2_H
#define CPUGREENSFUNCG2_H
//
// CPU Green's function for spreading charge on grid
// (c) Antti-Pekka Hynninen, March 2015
// aphynninen@hotmail.com
//
#include "CpuGrid.h"
#include "CpuGreensFunc.h"

template <typename T>
class CpuGreensFuncG2 : public CpuGreensFunc<T> {
private:
public:
  CpuGreensFuncG2() {}
  ~CpuGreensFuncG2() {}
  void mul(const int nfftX, const int nfftY, const int nfftZ,
	   const double boxx, const double boxy, const double boxz,
	   const double sigma, CpuGrid<T>& rho);
};
#endif //CPUGREENSFUNCG2_H
