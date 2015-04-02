
#include <cassert>
#include <cmath>
#include "CpuGSEr.h"

//
// Class creator
//
template <typename AT, typename CT>
CpuGSEr<AT, CT>::CpuGSEr(const int sizeX, const int sizeY, const int sizeZ,
			 const double sigma, const double kappa,
			 const double lambdaSigma, const double lambdaSigma1,
			 const double boxx, const double boxy, const double boxz,
			 const double ccelec) :
  CpuGSE<AT,CT>(sizeX, sizeY, sizeZ, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz, ccelec),
  rho(sizeX, sizeY, sizeZ), multiGridSolver(sizeX, sizeY, sizeZ) {
}

//
// Class destructor
//
template <typename AT, typename CT>
CpuGSEr<AT, CT>::~CpuGSEr() {
}

//
// First charge spreading
//
template <typename AT, typename CT>
void CpuGSEr<AT, CT>::spreadCharge1(const int numCoord, const xyzq_t<CT> *xyzq) {
  this->gaussCharge.spreadChargeToGrid(this->sigma1, this->lambdaSigma1*sqrt(2.0)*this->sigma1, numCoord, xyzq,
				       this->boxx, this->boxy, this->boxz, rho);
}

//
// Second charge spreading
//
template <typename AT, typename CT>
void CpuGSEr<AT, CT>::spreadCharge2() {
  this->gaussCharge.spreadChargeOnGrid(this->sigma2, this->lambdaSigma1*sqrt(2.0)*this->sigma2,
				       this->boxx, this->boxy, this->boxz, rho, this->rhoS);
}

//
// Solve Poisson equation
//
template <typename AT, typename CT>
void CpuGSEr<AT, CT>::solvePoisson() {
  this->phiM.clear();
  CpuGrid<CT> rhoTmp(this->rhoS.getSizeX(), this->rhoS.getSizeY(), this->rhoS.getSizeZ());
  rhoTmp.copy(this->rhoS);
  rhoTmp.scale(-4.0*pi);
  multiGridSolver.run(this->phiM, rhoTmp, this->boxx, this->boxy, this->boxz);
}

//
// Calculate energy
//
template <typename AT, typename CT>
AT CpuGSEr<AT, CT>::calculateEnergy() {
  return CpuGSE<AT,CT>::calculateEnergy(rho);
}

//
// Explicit instances of CpuGSEr
//
template class CpuGSEr<double, float>;
template class CpuGSEr<double, double>;
