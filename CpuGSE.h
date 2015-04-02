#ifndef CPUGSE_H
#define CPUGSE_H
//
// CPU version of Gaussian Split Ewald method
// This is the base class for real-space (CpuGSEr) and k-space (CpuGSEk) versions of GSE
// (c) Antti-Pekka Hynninen, Feb 2015
// aphynninen@hotmail.com
//
#include "TypesGSE.h"
#include "CpuGrid.h"
#include "CpuMultiGridSolver.h"
#include "CpuGaussCharge.h"

//
// AT  = Accumulation Type
// CT  = Calculation Type
//
template <typename AT, typename CT>
class CpuGSE {
protected:

  // Unit conversion, default is CHARMM value 332.0716
  const double ccelec;

  // Size of the simulation box in Angstroms
  double boxx;
  double boxy;
  double boxz;
  
  // Parameters
  double sigma;
  double kappa;
  double lambdaSigma;
  double lambdaSigma1;
  double sigma1;
  double sigma2;
  
  // Charge grid
  CpuGrid<CT> rhoS;

  // Grid-based electrostatic potential
  CpuGrid<CT> phiM;

  // Charge spreader
  CpuGaussCharge<AT,CT> gaussCharge;

  void setParam(const double sigma, const double kappa,
		const double lambdaSigma, const double lambdaSigma1);
  AT calculateEnergy(const CpuGrid<CT>& rhoIn);  
  
public:
  CpuGSE(const int sizeX, const int sizeY, const int sizeZ,
	 const double sigma, const double kappa,
	 const double lambdaSigma, const double lambdaSigma1,
	 const double boxx, const double boxy, const double boxz, const double ccelec=332.0716);
  ~CpuGSE();
  void printInfo();
  void setBoxSize(const double boxx, const double boxy, const double boxz);
  virtual void spreadCharge1(const int numCoord, const xyzq_t<CT> *xyzq)=0;
  virtual void solvePoisson()=0;
  void interpolateForce(const int numCoord, const xyzq_t<CT> *xyzq, AT *forceX, AT *forceY, AT *forceZ);
  virtual AT calculateEnergy()=0;
  CpuGrid<CT>& getRhoS() {return rhoS;}
  CpuGrid<CT>& getPhi() {return phiM;}
  double getSigma1() {return sigma1;}
  double getSigma2() {return sigma2;}
  double getLambdaSigma1() {return lambdaSigma1;}
};

#endif // CPUGSER_H
