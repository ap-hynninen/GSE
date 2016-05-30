#ifndef CPUEWALDDIRECT_H
#define CPUEWALDDIRECT_H
//
// Direct part of the Ewald summation
// Slow code, used as a reference only!
// (c) Antti-Pekka Hynninen May 2015, aphynninen@hotmail.com
//
#include "TypesGSE.h"

class CpuEwaldDirect {
private:

  // Unit conversion, default is CHARMM value 332.0716
  const double ccelec;
  
  // Width of the Gaussian
  double sigma;

  // Cut-off radius
  double rcut;
  
  // Box size
  double boxx;
  double boxy;
  double boxz;

  // Exclusion list
  int n_iblo14;
  int* iblo14;
  int n_inb14;
  int* inb14;
  
public:
  CpuEwaldDirect(const double sigma, const double rcut,
		 const double boxx, const double boxy, const double boxz,
		 const double ccelec=332.0716);
  ~CpuEwaldDirect();

  void setupExclusions(const int n_iblo14_in, const int* iblo14_in, const int* inb14_in);
  
  double calcForceEnergy(const int n, const xyzq_t<double>* xyzq,
			 double* fx, double* fy, double* fz);

};

#endif // CPUEWALDDIRECT_H
