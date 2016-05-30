#include <iostream>
#include <cassert>
#include <cmath>
#include "CpuEwaldDirect.h"

//
// Class constructor
//
CpuEwaldDirect::CpuEwaldDirect(const double sigma, const double rcut,
			       const double boxx, const double boxy, const double boxz,
			       const double ccelec) :
  sigma(sigma), rcut(rcut), ccelec(ccelec), boxx(boxx), boxy(boxy), boxz(boxz)
{
  n_iblo14 = 0;
  iblo14 = NULL;
  n_inb14 = 0;
  inb14 = NULL;
}

//
// Class destructor
//
CpuEwaldDirect::~CpuEwaldDirect() {
  if (iblo14 != NULL) delete [] iblo14;
  if (inb14 != NULL) delete [] inb14;
}

//
// Setup exclusions
//
void CpuEwaldDirect::setupExclusions(const int n_iblo14_in, const int* iblo14_in, const int* inb14_in) {
  n_iblo14 = n_iblo14_in;
  iblo14 = new int[n_iblo14];
  for (int i=0;i < n_iblo14;i++) iblo14[i] = iblo14_in[i];

  n_inb14 = iblo14[n_iblo14-1];
  inb14 = new int[n_inb14];
  // Copy and convert to 0-starting index
  for (int i=0;i < n_inb14;i++) inb14[i] = abs(inb14_in[i]) - 1;
}

//
// Calculates forces and energy (including self energy)
//
double CpuEwaldDirect::calcForceEnergy(const int n, const xyzq_t<double>* xyzq,
				       double* fx, double* fy, double* fz) {
  if (n_iblo14 > 0) assert(n_iblo14 == n);

  const double rcut2 = rcut*rcut;
  const double hboxx = 0.5*boxx;
  const double hboxy = 0.5*boxy;
  const double hboxz = 0.5*boxz;
  const double inv_sqrt2sigma = 1.0/(sqrt(2.0)*sigma);
  const double inv_2sigmasq = 1.0/(2.0*sigma*sigma);
  const double two_sqrtpi = 2.0/sqrt(pi_dbl);

  // Clear forces first
  for (int i=0;i < n;i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
  }

  double E = 0.0;
  for (int j=0;j < n-1;j++) {
    double xj = xyzq[j].x;
    double yj = xyzq[j].y;
    double zj = xyzq[j].z;
    double qj = xyzq[j].q*ccelec;
    int excl_start, excl_end;
    if (n_iblo14 > 0) {
      excl_start = (j >= 1) ? iblo14[j-1] : 0;
      excl_end = iblo14[j];
    }
    for (int i=j+1;i < n;i++) {
      bool exclude = false;
      if (n_iblo14 > 0) {
	for (int k=excl_start;k < excl_end;k++) {
	  if (i == inb14[k]) {
	    exclude = true;
	    break;
	  }
	}
      }
      if (exclude) continue;
      double xi = xyzq[i].x;
      double yi = xyzq[i].y;
      double zi = xyzq[i].z;
      double qi = xyzq[i].q;
      double dx = xi - xj;
      double dy = yi - yj;
      double dz = zi - zj;
      while (dx < -hboxx) dx += boxx;
      while (dx >  hboxx) dx -= boxx;
      while (dy < -hboxy) dy += boxy;
      while (dy >  hboxy) dy -= boxy;
      while (dz < -hboxz) dz += boxz;
      while (dz >  hboxz) dz -= boxz;
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < rcut2) {
	double qq = qj*qi;
	
	double r = sqrt(r2);
	double erfcval = erfc(r*inv_sqrt2sigma);
	double expval = exp(-r2*inv_2sigmasq);
	
	double f = -qq*(two_sqrtpi*inv_sqrt2sigma*expval + erfcval/r)/r2;
	double fijx = f*dx;
	double fijy = f*dy;
	double fijz = f*dz;

	fx[i] += fijx;
	fy[i] += fijy;
	fz[i] += fijz;
	fx[j] -= fijx;
	fy[j] -= fijy;
	fz[j] -= fijz;
	
	E += qq*erfcval/r;
      }
    }
  }

  double Eself = 0.0;
  for (int i=0;i < n;i++) {
    Eself += xyzq[i].q*xyzq[i].q;
  }
  Eself *= -ccelec/sqrt(2.0*pi_dbl*sigma*sigma);
  
  return (E + Eself);
}
