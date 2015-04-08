#ifndef CPUGAUSSCHARGE_H
#define CPUGAUSSCHARGE_H
//
// Class that performs Gaussian charge spreading and force interpolation
// (c) Antti-Pekka Hynninen, March 2015
//
#include <cstdio>
#include "TypesGSE.h"
#include "CpuGrid.h"

template <typename AT, typename CT>
class CpuGaussCharge {
private:
  // Size of the grid
  const int sizeX;
  const int sizeY;
  const int sizeZ;

  // Extra charge grids
  // NOTE: these are only allocated and used when method spreadChargeOnGrid is called
  CpuGrid<CT>* rhoX;
  CpuGrid<CT>* rhoXY;

  void interpolateLoop(const bool calcForce,
		       const double sigma, const double rcut,
		       const int numCoord, const xyzq_t<CT> *xyzq,
		       const double boxx, const double boxy, const double boxz,
		       const CpuGrid<CT>* phi,
		       const CpuGrid<CT>* Ex, const CpuGrid<CT>* Ey, const CpuGrid<CT>* Ez,
		       AT *resX, AT *resY, AT *resZ);

public:
  CpuGaussCharge(const int sizeX, const int sizeY, const int sizeZ);
  ~CpuGaussCharge();
  int calcSupportSize(const double rcut,
		      const double boxx, const double boxy, const double boxz);
  void calcElectricFieldOnGrid(const double boxx, const double boxy, const double boxz,
			       const CpuGrid<CT>& phi, CpuGrid<CT>& Ex, CpuGrid<CT>& Ey, CpuGrid<CT>& Ez);
  void interpolateForce(const double sigma, const double rcut,
			const int numCoord, const xyzq_t<CT> *xyzq,
			const double boxx, const double boxy, const double boxz,
			const CpuGrid<CT>& phi, AT *forceX, AT *forceY, AT *forceZ);
  void interpolateForce(const double sigma, const double rcut,
			const int numCoord, const xyzq_t<CT> *xyzq,
			const double boxx, const double boxy, const double boxz,
			const CpuGrid<CT>& Ex, const CpuGrid<CT>& Ey, const CpuGrid<CT>& Ez,
			AT *forceX, AT *forceY, AT *forceZ);
  void interpolateElectricField(const double sigma, const double rcut,
				const int numCoord, const xyzq_t<CT> *xyzq,
				const double boxx, const double boxy, const double boxz,
				const CpuGrid<CT>& ExM, const CpuGrid<CT>& EyM, const CpuGrid<CT>& EzM,
				AT* Ex, AT* Ey, AT* Ez);
  void interpolateElectricField(const double sigma, const double rcut,
				const int numCoord, const xyzq_t<CT> *xyzq,
				const double boxx, const double boxy, const double boxz,
				const CpuGrid<CT>& phi, AT* Ex, AT* Ey, AT* Ez);
  void calcDipoleSum(const double boxx, const double boxy, const double boxz,
		     const CpuGrid<CT>& rho, AT& dipSumX, AT& dipSumY, AT& dipSumZ);
  void spreadChargeToGrid(const double sigma, const double rcut,
			  const int numCoord, const xyzq_t<CT> *xyzq,
			  const double boxx, const double boxy, const double boxz,
			  CpuGrid<CT>& rho);
  void spreadChargeOnGrid(const double sigma, const double rcut,
			  const double boxx, const double boxy, const double boxz,
			  const CpuGrid<CT>& rho, CpuGrid<CT>& rhoS);
};
#endif //CPUGAUSSCHARGE_H
