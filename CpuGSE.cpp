
#include <cassert>
#include <cmath>
#include "CpuGSE.h"
#include "CpuGaussCharge.h"

//
// Class creator
//
template <typename AT, typename CT>
CpuGSE<AT, CT>::CpuGSE(const int sizeX, const int sizeY, const int sizeZ,
		       const double sigma, const double kappa,
		       const double lambdaSigma, const double lambdaSigma1,
		       const double boxx, const double boxy, const double boxz,
		       const double ccelec) :
  ccelec(ccelec), rhoS(sizeX, sizeY, sizeZ), phiM(sizeX, sizeY, sizeZ), gaussCharge(sizeX, sizeY, sizeZ) {

  setParam(sigma, kappa, lambdaSigma, lambdaSigma1);
  setBoxSize(boxx, boxy, boxz);
}

//
// Class destructor
//
template <typename AT, typename CT>
CpuGSE<AT, CT>::~CpuGSE() {
}

//
// Set parameters
//
template <typename AT, typename CT>
void CpuGSE<AT, CT>::setParam(const double sigma, const double kappa,
			      const double lambdaSigma, const double lambdaSigma1) {
  assert(kappa >= 0.0);
  assert(kappa <= 0.5);
  
  this->sigma = sigma;
  this->kappa = kappa;
  this->lambdaSigma = lambdaSigma;
  this->lambdaSigma1 = lambdaSigma1;
  sigma1 = sigma*sqrt(kappa);
  sigma2 = sigma*sqrt(1.0 - 2.0*kappa);
}

//
// Print info
//
template <typename AT, typename CT>
void CpuGSE<AT, CT>::printInfo() {
  std::cout << "sizeX=" << rhoS.getSizeX() << " sizeY=" << rhoS.getSizeY() << " sizeZ=" << rhoS.getSizeZ() << std::endl;
  std::cout << "sigma=" << sigma << " sigma1=" << sigma1 << " sigma2=" << sigma2 << " kappa=" << kappa << std::endl;
  std::cout << "lambdaSigma1=" << lambdaSigma1 << " lambdaSigma=" << lambdaSigma << std::endl;
}

//
// Set box size
//
template <typename AT, typename CT>
void CpuGSE<AT, CT>::setBoxSize(const double boxx, const double boxy, const double boxz) {
  this->boxx = boxx;
  this->boxy = boxy;
  this->boxz = boxz;
}

//
// Interpolate force from mesh-based potential
//
template <typename AT, typename CT>
void CpuGSE<AT, CT>::interpolateForce(const int numCoord, const xyzq_t<CT> *xyzq,
				      AT *forceX, AT *forceY, AT *forceZ) {
  this->gaussCharge.interpolateForce(this->sigma1, this->lambdaSigma1*sqrt(2.0)*this->sigma1,
				     numCoord, xyzq, this->boxx, this->boxy, this->boxz,
				     phiM, forceX, forceY, forceZ);
  for (int i=0;i < numCoord;i++) {
    forceX[i] *= ccelec;
    forceY[i] *= ccelec;
    forceZ[i] *= ccelec;
  }
}

//
// Calculate energy on mesh
//
template <typename AT, typename CT>
AT CpuGSE<AT, CT>::calculateEnergy(const CpuGrid<CT>& rhoIn) {
  // Sanity checks
  assert(rhoIn.getSizeX() == phiM.getSizeX());
  assert(rhoIn.getSizeY() == phiM.getSizeY());
  assert(rhoIn.getSizeZ() == phiM.getSizeZ());
  
  const int size  = rhoIn.getSize();
  const CT* rhoInData = rhoIn.getDataPointer();
  const CT* phiMdata = phiM.getDataPointer();
  
  AT energy = (AT)0;
  for (int i=0;i < size;i++) energy += (AT)(rhoInData[i]*phiMdata[i]);

  const AT hx = (AT)(boxx/(double)rhoIn.getSizeX());
  const AT hy = (AT)(boxy/(double)rhoIn.getSizeY());
  const AT hz = (AT)(boxz/(double)rhoIn.getSizeZ());
  energy *= ccelec*hx*hy*hz*(AT)0.5;
  
  return energy;
}

//
// Explicit instances of CpuGSE
//
template class CpuGSE<double, float>;
template class CpuGSE<double, double>;
