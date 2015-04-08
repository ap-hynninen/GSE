#include <iostream>
#include <cmath>
#include "CpuEwaldRecip.h"

//
// Class constructor
//
CpuEwaldRecip::CpuEwaldRecip(const double sigma, const int kcut,
			     const double boxx_in, const double boxy_in, const double boxz_in,
			     const double ccelec) :
  sigma(sigma), kcut(kcut), ccelec(ccelec) {
  setBoxsize(boxx_in, boxy_in, boxz_in);
}

//
// Class destructor
//
CpuEwaldRecip::~CpuEwaldRecip() {
}

//
// Calculates electric potential on a grid
//
void CpuEwaldRecip::calcElectricPotential(const int n, const xyzq_t<double>* xyzq, CpuGrid<double>& phi) {
  const int sizeX = phi.getSizeX();
  const int sizeY = phi.getSizeY();
  const int sizeZ = phi.getSizeZ();
  const int size = sizeX*sizeY*sizeZ;
  double* phiData = phi.getDataPointer();
  const double hx = boxx/(double)sizeX;
  const double hy = boxy/(double)sizeY;
  const double hz = boxz/(double)sizeZ;

  const double inv_boxx = 1.0/boxx;
  const double inv_boxy = 1.0/boxy;
  const double inv_boxz = 1.0/boxz;
  const double one_4sigma2 = 1.0/(4.0*sigma*sigma);
  const double pref = 2.0*M_PI*inv_boxx*inv_boxy*inv_boxz*ccelec;

  phi.clear();
  
  int kcut2 = kcut*kcut;
  double E = 0.0;
  for (int kz=-kcut;kz <= kcut;kz++) {
    double kzv = 2.0*M_PI*((double)kz)*inv_boxz;
    for (int ky=-kcut;ky <= kcut;ky++) {
      double kyv = 2.0*M_PI*((double)ky)*inv_boxy;
      for (int kx=-kcut;kx <= kcut;kx++) {
	int k2 = kx*kx + ky*ky + kz*kz;
	if ((kx == 0 && ky == 0 && kz == 0) || k2 > kcut2) continue;
	double kxv = 2.0*M_PI*((double)kx)*inv_boxx;
	double ksq = kxv*kxv + kyv*kyv + kzv*kzv;
	double Sk_re = 0.0;
	double Sk_im = 0.0;
	for (int i=0;i < n;i++) {
	  double xi = xyzq[i].x;
	  double yi = xyzq[i].y;
	  double zi = xyzq[i].z;
	  double kr = kxv*xi + kyv*yi + kzv*zi;
	  double qi = xyzq[i].q;
	  double cos_kr = cos(kr);
	  double sin_kr = sin(kr);
	  Sk_re += qi*cos_kr;
	  Sk_im += qi*sin_kr;
	}
	double pref_exp_k = pref*exp(-ksq*one_4sigma2)/ksq;
	int p = 0;
	for (int iz=0;iz < sizeZ;iz++) {
	  for (int iy=0;iy < sizeY;iy++) {
	    for (int ix=0;ix < sizeX;ix++,p++) {
	      double kr = kxv*hx*(double)ix + kyv*hy*(double)iy + kzv*hz*(double)iz;
	      double cos_kr = cos(kr);
	      double sin_kr = sin(kr);
	      phiData[p] += cos_kr*Sk_re - sin_kr*Sk_im;
	    }
	  }
	}
      }
    }
  }
  
}

//
// Calculates forces and energy. Nk = number of k vectors to include in each direction
//
double CpuEwaldRecip::calcForceEnergy(const int n, const xyzq_t<double>* xyzq,
				      double* fx, double* fy, double* fz) {

  double inv_boxx = 1.0/boxx;
  double inv_boxy = 1.0/boxy;
  double inv_boxz = 1.0/boxz;
  double half_sigmasq = 0.5*sigma*sigma;
  double pref = 2.0*M_PI*inv_boxx*inv_boxy*inv_boxz*ccelec;

  // Clear forces first
  for (int i=0;i < n;i++) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
  }
  
  int kcut2 = kcut*kcut;
  double E = 0.0;
  for (int kz=-kcut;kz <= kcut;kz++) {
    double kzv = 2.0*M_PI*((double)kz)*inv_boxz;
    for (int ky=-kcut;ky <= kcut;ky++) {
      double kyv = 2.0*M_PI*((double)ky)*inv_boxy;
      for (int kx=-kcut;kx <= kcut;kx++) {
	int k2 = kx*kx + ky*ky + kz*kz;
	if ((kx == 0 && ky == 0 && kz == 0) || k2 > kcut2) continue;
	double kxv = 2.0*M_PI*((double)kx)*inv_boxx;
	double ksq = kxv*kxv + kyv*kyv + kzv*kzv;
	double rhok_re = 0.0;
	double rhok_im = 0.0;
	for (int i=0;i < n;i++) {
	  double xi = xyzq[i].x;
	  double yi = xyzq[i].y;
	  double zi = xyzq[i].z;
	  double kr = kxv*xi + kyv*yi + kzv*zi;
	  double qi = xyzq[i].q;
	  double sin_kr = sin(kr);
	  double cos_kr = cos(kr);
	  rhok_re += qi*cos_kr;
	  rhok_im += qi*sin_kr;
	}
	double pref_exp_k = pref*exp(-ksq*half_sigmasq)/ksq;
	for (int i=0;i < n;i++) {
	  double xi = xyzq[i].x;
	  double yi = xyzq[i].y;
	  double zi = xyzq[i].z;
	  double kr = kxv*xi + kyv*yi + kzv*zi;
	  double qi = xyzq[i].q;
	  double sin_kr = sin(kr);
	  double cos_kr = cos(kr);
	  double fac = -pref_exp_k*qi*2.0;
	  fx[i] += fac*(rhok_re*(-kxv*sin_kr) + rhok_im*kxv*cos_kr);
	  fy[i] += fac*(rhok_re*(-kyv*sin_kr) + rhok_im*kyv*cos_kr);
	  fz[i] += fac*(rhok_re*(-kzv*sin_kr) + rhok_im*kzv*cos_kr);
	}
	E += pref_exp_k*(rhok_re*rhok_re + rhok_im*rhok_im);
      }
    }
  }
 
  return E;
}
