#ifndef CPUGSER_H
#define CPUGSER_H
//
// CPU version of Gaussian Split Ewald method, real-space version
// (c) Antti-Pekka Hynninen, Feb 2015
// aphynninen@hotmail.com
//
#include "TypesGSE.h"
#include "CpuGrid.h"
#include "CpuMultiGridSolver.h"
#include "CpuGSE.h"

//
// AT  = Accumulation Type
// CT  = Calculation Type
//
template <typename AT, typename CT>
class CpuGSEr : public CpuGSE<AT,CT> {
private:
    
  // Multi grid solver for Poisson equation
  CpuMultiGridSolver<CT> multiGridSolver;

  // Initial charge density
  CpuGrid<CT> rho;

public:
  CpuGSEr(const int sizeX, const int sizeY, const int sizeZ,
	 const double sigma, const double kappa,
	 const double lambdaSigma, const double lambdaSigma1,
	 const double boxx, const double boxy, const double boxz, const double ccelec=332.0716);
  ~CpuGSEr();
  void spreadCharge1(const int numCoord, const xyzq_t<CT> *xyzq);
  void spreadCharge2();
  void solvePoisson();
  AT calculateEnergy();
};

#endif // CPUGSER_H
