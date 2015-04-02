#ifndef CPUGREENSFUNCGSE_H
#define CPUGREENSFUNCGSE_H
//
// CPU Green's function for k-space GSE.
// (c) Antti-Pekka Hynninen, March 2015
// aphynninen@hotmail.com
//
#include "CpuGrid.h"
#include "CpuGreensFunc.h"

template <typename T>
class CpuGreensFuncGSE : public CpuGreensFunc<T> {
private:
public:
  CpuGreensFuncGSE() {}
  ~CpuGreensFuncGSE() {}
  void mul(const int nfftX, const int nfftY, const int nfftZ,
	   const double boxx, const double boxy, const double boxz,
	   const double sigma, CpuGrid<T>& rho);
};
#endif //CPUGREENSFUNCGSE_H
