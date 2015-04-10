#include <cassert>
#include "TypesGSE.h"
#include "CpuGreensFuncG2.h"

//
// Multiply Fourier transformed charge density by the Green's function
//
template <typename T>
void CpuGreensFuncG2<T>::mul(const int nfftX, const int nfftY, const int nfftZ,
			     const double boxx, const double boxy, const double boxz,
			     const double sigma, CpuGrid<T>& rho) {
  // Sanity checks
  assert((nfftX/2+1)*2 == rho.getSizeX());
  assert(nfftY == rho.getSizeY());
  assert(nfftZ == rho.getSizeZ());

  const int sizeX = nfftX/2+1;
  const int sizeY = nfftY;
  const int sizeZ = nfftZ;
  
  const int halfNfftX = nfftX/2 + (nfftX % 2);
  const int halfNfftY = nfftY/2 + (nfftY % 2);
  const int halfNfftZ = nfftZ/2 + (nfftZ % 2);

  const T twopi_boxx = (T)(2.0*pi_dbl/boxx);
  const T twopi_boxy = (T)(2.0*pi_dbl/boxy);
  const T twopi_boxz = (T)(2.0*pi_dbl/boxz);
  //const T inv_vol = (T)(1.0/(boxx*boxy*boxz));
  const T inv_vol = (T)(1.0/(double)(nfftX*nfftY*nfftZ));
  const T sigmafac = (T)(sigma*sigma/2.0);

  T* rhoData = rho.getDataPointer();
  int p = 0;
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++,p+=2) {
	//if (p == 0) {
	//  rhoData[0] = (T)0.0;
	//  rhoData[1] = (T)0.0;
	//  continue;
	//}
	//int nx = (ix >= halfNfftX) ? ix-nfftX : ix;
	//int ny = (iy >= halfNfftY) ? iy-nfftY : iy;
	//int nz = (iz >= halfNfftZ) ? iz-nfftZ : iz;

	int kx = (ix >= halfNfftX) ? ix-nfftX : ix;
	int ky = (iy >= halfNfftY) ? iy-nfftY : iy;
	int kz = (iz >= halfNfftZ) ? iz-nfftZ : iz;
	
	T kxd = ((T)kx)*twopi_boxx;
	T kyd = ((T)ky)*twopi_boxy;
	T kzd = ((T)kz)*twopi_boxz;
	T kd2 = kxd*kxd + kyd*kyd + kzd*kzd;
	T fac = inv_vol*exp(-sigmafac*kd2);
	rhoData[p]   *= fac;
	rhoData[p+1] *= fac;
      }
    }
  }

}

//
// Explicit instances of CpuGreensFuncG2
//
template class CpuGreensFuncG2<float>;
template class CpuGreensFuncG2<double>;

