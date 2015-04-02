#ifndef CPUGSEK_H
#define CPUGSEK_H
//
// CPU version of Gaussian Split Ewald method, k-space version
// (c) Antti-Pekka Hynninen, Feb 2015
// aphynninen@hotmail.com
//
#include "TypesGSE.h"
#include "CpuGrid.h"
#include "CpuGSE.h"
#include "CpuFFTSolver.h"
#include "CpuGreensFuncGSE.h"

//
// AT  = Accumulation Type
// CT  = Calculation Type
//
template <typename AT, typename CT>
class CpuGSEk : public CpuGSE<AT,CT> {
private:
  
  // FFT solver for Poisson equation
  CpuFFTSolver<CT> fftSolver;

  CpuGreensFuncGSE<CT> greensFunc;
  
public:
  CpuGSEk(const int sizeX, const int sizeY, const int sizeZ,
	 const double sigma, const double kappa,
	 const double lambdaSigma, const double lambdaSigma1,
	 const double boxx, const double boxy, const double boxz, const double ccelec=332.0716);
  ~CpuGSEk();
  void spreadCharge1(const int numCoord, const xyzq_t<CT> *xyzq);
  void spreadCharge2(CpuGrid<CT>& rhoOut);
  void solvePoisson();
  AT calculateEnergy();
};

#endif // CPUGSEK_H
