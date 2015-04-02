#ifndef CPUGSEU_H
#define CPUGSEU_H
//
// CPU version of Gaussian Split Ewald method, u-series version
// (c) Antti-Pekka Hynninen, March 2015
// aphynninen@hotmail.com
//
#include <vector>
#include "TypesGSE.h"
#include "CpuGrid.h"
#include "CpuGSE.h"

//
// AT  = Accumulation Type
// CT  = Calculation Type
//
template <typename AT, typename CT>
class CpuGSEu : public CpuGSE<AT,CT> {
private:

  // Number of terms in the u-series
  const int N;

  // constant for u-series
  const double b;

  // Set of intermediate charge densities, size N
  std::vector<CpuGrid<CT>*> rhoG;
  
public:
  CpuGSEu(const int N, const double u,
	  const int sizeX, const int sizeY, const int sizeZ,
	  const double sigma, const double kappa,
	  const double lambdaSigma, const double lambdaSigma1,
	  const double boxx, const double boxy, const double boxz, const double ccelec=332.0716);
  ~CpuGSEu();
  // Overriding
  void interpolateForce(const int numCoord, const xyzq_t<CT> *xyzq, AT *forceX, AT *forceY, AT *forceZ);
  void spreadCharge1(const int numCoord, const xyzq_t<CT> *xyzq);
  void solvePoisson();
  AT calculateEnergy();
};

#endif // CPUGSEU_H
