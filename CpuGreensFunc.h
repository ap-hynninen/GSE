#ifndef CPUGREENSFUNC_H
#define CPUGREENSFUNC_H
//
// CPU Green's function base class. This class multiplies k-space charge density by the Green's function
// (c) Antti-Pekka Hynninen, March 2015
// aphynninen@hotmail.com
//
#include "CpuGrid.h"

template <typename T>
class CpuGreensFunc {
private:
public:
  CpuGreensFunc() {}
  ~CpuGreensFunc() {}
  virtual void mul(const int nfftX, const int nfftY, const int nfftZ,
		   const double boxx, const double boxy, const double boxz,
		   const double sigma, CpuGrid<T>& rho)=0;
};
#endif //CPUGREENSFUNC_H
