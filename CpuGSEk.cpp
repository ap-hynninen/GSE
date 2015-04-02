
#include <cassert>
#include <cmath>
#include "CpuGSEk.h"

//
// Class creator
//
template <typename AT, typename CT>
CpuGSEk<AT, CT>::CpuGSEk(const int sizeX, const int sizeY, const int sizeZ,
			 const double sigma, const double kappa,
			 const double lambdaSigma, const double lambdaSigma1,
			 const double boxx, const double boxy, const double boxz,
			 const double ccelec) :
  CpuGSE<AT,CT>(sizeX, sizeY, sizeZ, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz, ccelec),
  fftSolver(greensFunc, sizeX, sizeY, sizeZ) {
}

//
// Class destructor
//
template <typename AT, typename CT>
CpuGSEk<AT, CT>::~CpuGSEk() {
}

//
// First charge spreading
//
template <typename AT, typename CT>
void CpuGSEk<AT, CT>::spreadCharge1(const int numCoord, const xyzq_t<CT> *xyzq) {
  this->gaussCharge.spreadChargeToGrid(this->sigma1, this->lambdaSigma1*sqrt(2.0)*this->sigma1, numCoord, xyzq,
				       this->boxx, this->boxy, this->boxz, this->rhoS);
}

//
// Second charge spreading
//
template <typename AT, typename CT>
void CpuGSEk<AT, CT>::spreadCharge2(CpuGrid<CT>& rhoOut) {
  this->gaussCharge.spreadChargeOnGrid(this->sigma2, this->lambdaSigma1*sqrt(2.0)*this->sigma2,
				       this->boxx, this->boxy, this->boxz, this->rhoS, rhoOut);
}

//
// Solve Poisson equation
//
template <typename AT, typename CT>
void CpuGSEk<AT, CT>::solvePoisson() {
  fftSolver.run(this->phiM, this->rhoS, this->sigma2, this->boxx, this->boxy, this->boxz);
}

//
// Calculate energy
//
template <typename AT, typename CT>
AT CpuGSEk<AT, CT>::calculateEnergy() {
  return CpuGSE<AT,CT>::calculateEnergy(this->rhoS);
}

//
// Explicit instances of CpuGSEk
//
template class CpuGSEk<double, float>;
template class CpuGSEk<double, double>;
