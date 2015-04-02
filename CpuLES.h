#ifndef CPULES_H
#define CPULES_H
//
// Local electrostatics implementation on CPU
// 
// Based on the algorithm by Jorg Rottler et al:
// Jorg Rottler, JCP 127, 134104 (2007)
//
// (c) Antti-Pekka Hynninen 2014
//
#include "CpuGrid.h"
#include "CpuGaussCharge.h"

template<typename AT, typename CT>
class CpuLES {
 private:

  enum {DX, DY, DZ};
  
  // Unit conversion, default is CHARMM value 332.0716
  const double ccelec;

  // Width of the Gaussian
  double sigma;
  
  // Size of the simulation box in Angstroms
  double boxx;
  double boxy;
  double boxz;

  // Electric field components
  CpuGrid<CT> Ex;
  CpuGrid<CT> Ey;
  CpuGrid<CT> Ez;

  // Magnetic field components
  CpuGrid<CT> Bx;
  CpuGrid<CT> By;
  CpuGrid<CT> Bz;

  // Charge density
  CpuGrid<CT> rho;
  CpuGrid<CT> rhoTmp;

  // Charge spreader
  CpuGaussCharge<AT,CT> gaussCharge;

  void calcDipoleSum(const int numCoord, const xyzq_t<CT> *xyzq,
		     AT& dipSumX, AT& dipSumY, AT& dipSumZ);

 public:
  CpuLES(const int sizeX, const int sizeY, const int sizeZ, const double sigma,
	 const double boxx, const double boxy, const double boxz,
	 const double ccelec=332.0716);
  ~CpuLES();

  double calcTotalEnergy();
  void interpolateForceVW(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq,
			  AT *forceX, AT *forceY, AT *forceZ);
  void interpolateForceEF(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq,
			  AT *forceX, AT *forceY, AT *forceZ);
  void interpolateElectricField(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq,
				AT *ExPart, AT *EyPart, AT *EzPart);
  void spreadCharge1(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq);
  void spreadCharge2(const double sigma2, const double lambdaSigma1);
  void initElectricFieldHolm();
  void initElectricField();
  double checkGaussLaw();
  void integrate(const CT c, const CT dt, const CT gamma2);
  void clearMagneticField();
  void setElectricField(const CpuGrid<CT>& ExIn, const CpuGrid<CT>& EyIn, const CpuGrid<CT>& EzIn);

  CT curl(const int curlType, const int ix, const int iy, const int iz,
	  const int sizeX, const int sizeY, const int sizeZ,
	  const CT inv_h,
	  const CT* DxData, const CT* DyData, const CT* DzData,
	  CT& valx, CT& valy, CT& valz);
  CT maxCurl(const int curlType, const CpuGrid<CT>& Dx, const CpuGrid<CT>& Dy, const CpuGrid<CT>& Dz);
  CT maxCurlB();
  CT maxCurlE();
  double calcDipoleEnergy(const int numCoord, const xyzq_t<CT>* xyzq);
  
  CpuGrid<CT>& getEx() {return Ex;}
  CpuGrid<CT>& getEy() {return Ey;}
  CpuGrid<CT>& getEz() {return Ez;}
  CpuGrid<CT>& getRho() {return rho;}
};

#endif // CPULES_H
