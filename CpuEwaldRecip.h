#ifndef CPUEWALDRECIP_H
#define CPUEWALDRECIP_H
//
// Calculates the reciprocal part of the Ewald summation.
// Used mainly as a reference.
//
// (c) Antti-Pekka Hynninen 2014
//
#include "TypesGSE.h"
#include "CpuGrid.h"

class CpuEwaldRecip {
 private:
  // Unit conversion, default is CHARMM value 332.0716
  const double ccelec;
  
  // Width of the Gaussian
  double sigma;

  // Cut-off for k-vectors
  int kcut;
  
  // Box size
  double boxx;
  double boxy;
  double boxz;

 public:
  CpuEwaldRecip(const double sigma, const int kcut,
		const double boxx_in, const double boxy_in, const double boxz_in,
		const double ccelec=332.0716);
  ~CpuEwaldRecip();

  void calcElectricPotential(const int n, const xyzq_t<double>* xyzq, CpuGrid<double>& phi);
  double calcForceEnergy(const int n, const xyzq_t<double>* xyzq,
			 double* fx, double* fy, double* fz);

  void setSigma(const double sigma_in) {this->sigma = sigma_in;}
  
  void setBoxsize(const double boxx_in, const double boxy_in, const double boxz_in) {
    this->boxx = boxx_in;
    this->boxy = boxy_in;
    this->boxz = boxz_in;
  }
  
};

#endif // CPUEWALDRECIP_H
