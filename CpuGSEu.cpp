
#include <cassert>
#include <cmath>
#include "CpuGSEu.h"

//
// Class creator
//
template <typename AT, typename CT>
CpuGSEu<AT, CT>::CpuGSEu(const int N, const double b,
			 const int sizeX, const int sizeY, const int sizeZ,
			 const double sigma, const double kappa,
			 const double lambdaSigma, const double lambdaSigma1,
			 const double boxx, const double boxy, const double boxz,
			 const double ccelec) :
  N(N), b(b), CpuGSE<AT,CT>(sizeX, sizeY, sizeZ, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz, ccelec) {
  assert(N > 0);
  assert(b > 1.0);
  
  rhoG.reserve(N);
  for (int j=0;j < N;j++) {
    rhoG.push_back(new CpuGrid<CT>(sizeX, sizeY, sizeZ));
  }
}

//
// Class destructor
//
template <typename AT, typename CT>
CpuGSEu<AT, CT>::~CpuGSEu() {
  for (int j=0;j < N;j++) {
    delete rhoG.at(j);
  }
}

//
// First charge spreading
//
template <typename AT, typename CT>
void CpuGSEu<AT, CT>::spreadCharge1(const int numCoord, const xyzq_t<CT> *xyzq) {
  this->gaussCharge.spreadChargeToGrid(this->sigma1, this->lambdaSigma1*sqrt(2.0)*this->sigma1, numCoord, xyzq,
				       this->boxx, this->boxy, this->boxz, this->rhoS);
}

//
// Solve Poisson equation using u-series
//
template <typename AT, typename CT>
void CpuGSEu<AT, CT>::solvePoisson() {
  // Compute Convolutions
  for (int j=0;j < N;j++) {
    double sigmaj = sqrt(this->sigma*this->sigma*pow(b,2*j) - 2*this->sigma1*this->sigma1);
    //printf("sigmaj = %lf\n",sigmaj);
    this->gaussCharge.spreadChargeOnGrid(sigmaj, this->lambdaSigma1*sqrt(2.0)*sigmaj,
					 this->boxx, this->boxy, this->boxz, this->rhoS, *rhoG.at(j));
    double fac = pow(2*pi_dbl*sigmaj*sigmaj,3.0/2.0)/pow(sigmaj*sigmaj,1.0/2.0);
    rhoG.at(j)->scale(fac);
  }
  // Combine convolutions to a single grid
  this->phiM.clear();
  for (int j=N-1;j >=0 ;j--) {
    char filename[64];
    sprintf(filename,"rhoG%d.txt",j);
    rhoG.at(j)->save(filename);
    this->phiM.add(*rhoG.at(j));
  }  
}

//
// Interpolate force from mesh-based potential
//
template <typename AT, typename CT>
void CpuGSEu<AT, CT>::interpolateForce(const int numCoord, const xyzq_t<CT> *xyzq,
				       AT *forceX, AT *forceY, AT *forceZ) {
  this->gaussCharge.interpolateForce(this->sigma1, this->lambdaSigma1*sqrt(2.0)*this->sigma1,
				     numCoord, xyzq, this->boxx, this->boxy, this->boxz,
				     this->phiM, forceX, forceY, forceZ);
  AT pref = 2.0*log(b)*this->ccelec;
  for (int i=0;i < numCoord;i++) {
    forceX[i] *= pref;
    forceY[i] *= pref;
    forceZ[i] *= pref;
  }
}

//
// Calculate energy
//
template <typename AT, typename CT>
AT CpuGSEu<AT, CT>::calculateEnergy() {
  return 2.0*log(b)*CpuGSE<AT,CT>::calculateEnergy(this->rhoS);
}

//
// Explicit instances of CpuGSEu
//
template class CpuGSEu<double, float>;
template class CpuGSEu<double, double>;
