#ifndef CUDAGAUSSCHARGE_H
#define CUDAGAUSSCHARGE_H
//
// Class that performs Gaussian charge spreading and force interpolation on GPU
// (c) Antti-Pekka Hynninen, March 2015
//
#include <cstdio>
#include "TypesGSE.h"
#include "CudaGrid.h"

template <typename AT, typename CT>
class CudaGaussCharge {
private:
  // Size of the grid
  const int sizeX;
  const int sizeY;
  const int sizeZ;

public:
  CudaGaussCharge(const int sizeX, const int sizeY, const int sizeZ);
  ~CudaGaussCharge();
  void spreadChargeToGrid(const double sigma, const double rcut,
			  const int numCoord, const xyzq_t<CT> *xyzq,
			  const double boxx, const double boxy, const double boxz,
			  CudaGrid<CT>& rho);
};
#endif //CUDAGAUSSCHARGE_H
