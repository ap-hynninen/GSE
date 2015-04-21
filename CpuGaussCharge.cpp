#include <cassert>
#include <cmath>
#include "CpuGaussCharge.h"

//
// Class creator
//
template <typename AT, typename CT>
CpuGaussCharge<AT, CT>::CpuGaussCharge(const int sizeX, const int sizeY, const int sizeZ) :
  sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ) {
  rhoX = NULL;
  rhoXY = NULL;
}

//
// Class destructor
//
template <typename AT, typename CT>
CpuGaussCharge<AT, CT>::~CpuGaussCharge() {
  if (rhoX != NULL) delete rhoX;
  if (rhoXY != NULL) delete rhoXY;
}

//
// Calculate size of the support
//
template <typename AT, typename CT>
int CpuGaussCharge<AT, CT>::calcSupportSize(const double rcut,
					    const double boxx, const double boxy, const double boxz) {
  const double hxd = boxx/(double)sizeX;
  const double hyd = boxy/(double)sizeY;
  const double hzd = boxz/(double)sizeZ;
  const int nx = (int)ceil(rcut/hxd);
  const int ny = (int)ceil(rcut/hyd);
  const int nz = (int)ceil(rcut/hzd);

  const CT hx = (CT)hxd;
  const CT hy = (CT)hyd;
  const CT hz = (CT)hzd;
  const CT rcut2 = (CT)(rcut*rcut);

  int nsup = 0;
  for (int tz=-nz;tz <= nz;tz++) {
    CT dz = ((CT)(abs(tz)-1))*hz;
    for (int ty=-ny;ty <= ny;ty++) {
      CT dy = ((CT)(abs(ty)-1))*hy;
      for (int tx=-nx;tx <= nx;tx++) {
	CT dx = ((CT)(abs(tx)-1))*hx;
	CT r2 = dx*dx + dy*dy + dz*dz;
	if (r2 < rcut2) nsup++;
      }
    }
  }

  return nsup;
}

//
// Calculates electric field on grid
//
template <typename AT, typename CT>
void CpuGaussCharge<AT, CT>::calcElectricFieldOnGrid(const double boxx, const double boxy, const double boxz,
						     const CpuGrid<CT>& phi, CpuGrid<CT>& Ex, CpuGrid<CT>& Ey, CpuGrid<CT>& Ez) {
  // Sanity checks
  assert(sizeX == phi.getSizeX());
  assert(sizeY == phi.getSizeY());
  assert(sizeZ == phi.getSizeZ());
  assert(sizeX == Ex.getSizeX());
  assert(sizeY == Ex.getSizeY());
  assert(sizeZ == Ex.getSizeZ());
  assert(sizeX == Ey.getSizeX());
  assert(sizeY == Ey.getSizeY());
  assert(sizeZ == Ey.getSizeZ());
  assert(sizeX == Ez.getSizeX());
  assert(sizeY == Ez.getSizeY());
  assert(sizeZ == Ez.getSizeZ());
  const double hxd = boxx/(double)sizeX;
  const double hyd = boxy/(double)sizeY;
  const double hzd = boxz/(double)sizeZ;
  //const CT inv_2hx = (CT)(0.5/hxd);
  //const CT inv_2hy = (CT)(0.5/hyd);
  //const CT inv_2hz = (CT)(0.5/hzd);
  const CT inv_hx = (CT)(1.0/hxd);
  const CT inv_hy = (CT)(1.0/hyd);
  const CT inv_hz = (CT)(1.0/hzd);
  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int p = ix + (iy + iz*sizeY)*sizeX;
	ExData[p] = -(phi.getDataValue(ix+1, iy, iz) - phi.getDataValue(ix, iy, iz))*inv_hx;
	EyData[p] = -(phi.getDataValue(ix, iy+1, iz) - phi.getDataValue(ix, iy, iz))*inv_hy;
	EzData[p] = -(phi.getDataValue(ix, iy, iz+1) - phi.getDataValue(ix, iy, iz))*inv_hz;
      }
    }
  }
}

//
// Interpolation loop for force and electric field computation
//
template <typename AT, typename CT>
void CpuGaussCharge<AT, CT>::interpolateLoop(const double sigma, const double rcut,
					     const int numCoord, const xyzq_t<CT> *xyzq,
					     const double boxx, const double boxy, const double boxz,
					     const CpuGrid<CT>& phi, AT *resX, AT *resY, AT *resZ) {
  // Sanity checks
  assert(sizeX == phi.getSizeX());
  assert(sizeY == phi.getSizeY());
  assert(sizeZ == phi.getSizeZ());

  const double hxd = boxx/(double)sizeX;
  const double hyd = boxy/(double)sizeY;
  const double hzd = boxz/(double)sizeZ;
  const int nx = (int)ceil(rcut/hxd);
  const int ny = (int)ceil(rcut/hyd);
  const int nz = (int)ceil(rcut/hzd);
  const CT inv_hx = (CT)(1.0/hxd);
  const CT inv_hy = (CT)(1.0/hyd);
  const CT inv_hz = (CT)(1.0/hzd);
  const CT hx = (CT)hxd;
  const CT hy = (CT)hyd;
  const CT hz = (CT)hzd;
  const CT rcut2 = (CT)(rcut*rcut);
  const CT inv_2sigmasq = (CT)(1.0/(2.0*sigma*sigma));
  const CT pref = (CT)(pow(2.0*pi_dbl,-3.0/2.0)*hxd*hyd*hzd*pow(sigma,-5.0));
  const CT* phiData = phi.getDataPointer();

  printf("CpuGaussCharge::interpolateLoop, rcut=%lf nx=%d support=%d\n",rcut,
	 nx,calcSupportSize(rcut, boxx, boxy, boxz));

  for (int i=0;i < numCoord;i++) {
    CT x = xyzq[i].x;
    CT y = xyzq[i].y;
    CT z = xyzq[i].z;
    CT fac = xyzq[i].q*pref;
    int ix = (int)round(x*inv_hx);
    int iy = (int)round(y*inv_hy);
    int iz = (int)round(z*inv_hz);
    AT rx = (AT)0.0;
    AT ry = (AT)0.0;
    AT rz = (AT)0.0;
    for (int tz=iz-nz;tz <= iz+nz;tz++) {
      CT zp = ((CT)tz)*hz;
      CT dz = z-zp;
      int pz = (tz + sizeZ) % sizeZ;
      for (int ty=iy-ny;ty <= iy+ny;ty++) {
	CT yp = ((CT)ty)*hy;
	CT dy = y-yp;
	int py = (ty + sizeY) % sizeY;
	for (int tx=ix-nx;tx <= ix+nx;tx++) {
	  CT xp = ((CT)tx)*hx;
	  CT dx = x-xp;
	  int px = (tx + sizeX) % sizeX;
	  CT r2 = dx*dx + dy*dy + dz*dz;
	  if (r2 < rcut2) {
	    int p = px + (py + pz*sizeY)*sizeX;
	    CT r = fac*exp(-r2*inv_2sigmasq)*phiData[p];
	    rx += (AT)(r*dx);
	    ry += (AT)(r*dy);
	    rz += (AT)(r*dz);
	  }
	}
      }
    }
    resX[i] = rx;
    resY[i] = ry;
    resZ[i] = rz;
  }
}

//
// Interpolate force from grid
//
template <typename AT, typename CT>
void CpuGaussCharge<AT, CT>::interpolateForce(const double sigma, const double rcut,
					      const int numCoord, const xyzq_t<CT> *xyzq,
					      const double boxx, const double boxy, const double boxz,
					      const CpuGrid<CT>& phi, AT *forceX, AT *forceY, AT *forceZ) {
  interpolateLoop(sigma, rcut, numCoord, xyzq, boxx, boxy, boxz, phi, forceX, forceY, forceZ);
}

//
// Calculate dipole sum of the charge density
//
template <typename AT, typename CT>
void CpuGaussCharge<AT, CT>::calcDipoleSum(const double boxx, const double boxy, const double boxz,
					   const CpuGrid<CT>& rho, AT& dipSumX, AT& dipSumY, AT& dipSumZ) {
  // Sanity checks
  assert(sizeX == rho.getSizeX());
  assert(sizeY == rho.getSizeY());
  assert(sizeZ == rho.getSizeZ());
  
  const double hxd = boxx/(double)sizeX;
  const double hyd = boxy/(double)sizeY;
  const double hzd = boxz/(double)sizeZ;
  const CT inv_hx = (CT)(1.0/hxd);
  const CT inv_hy = (CT)(1.0/hyd);
  const CT inv_hz = (CT)(1.0/hzd);
  const CT hx = (CT)hxd;
  const CT hy = (CT)hyd;
  const CT hz = (CT)hzd;

  dipSumX = (AT)0.0;
  dipSumY = (AT)0.0;
  dipSumZ = (AT)0.0;
  
  const CT* rhodata = rho.getDataPointer();
  for (int iz=0;iz < sizeZ;iz++) {
    CT z = hz*(CT)iz;
    for (int iy=0;iy < sizeY;iy++) {
      CT y = hy*(CT)iy;
      for (int ix=0;ix < sizeX;ix++) {
	CT x = hx*(CT)ix;
	int p = ix + (iy + iz*sizeY)*sizeX;
	CT val = rhodata[p]*hx*hy*hz;
	dipSumX += (AT)(val*x);
	dipSumY += (AT)(val*y);
	dipSumZ += (AT)(val*z);
      }
    }
  }

}

//
// Particle to grid charge spreading
//
template <typename AT, typename CT>
void CpuGaussCharge<AT, CT>::spreadChargeToGrid(const double sigma, const double rcut,
						const int numCoord, const xyzq_t<CT> *xyzq,
						const double boxx, const double boxy, const double boxz,
						CpuGrid<CT>& rho) {
  // Sanity checks
  assert(sizeX == rho.getSizeX());
  assert(sizeY == rho.getSizeY());
  assert(sizeZ == rho.getSizeZ());
  
  const double hxd = boxx/(double)sizeX;
  const double hyd = boxy/(double)sizeY;
  const double hzd = boxz/(double)sizeZ;
  const int nx = (int)ceil(rcut/hxd);
  const int ny = (int)ceil(rcut/hyd);
  const int nz = (int)ceil(rcut/hzd);
  const CT inv_hx = (CT)(1.0/hxd);
  const CT inv_hy = (CT)(1.0/hyd);
  const CT inv_hz = (CT)(1.0/hzd);
  const CT hx = (CT)hxd;
  const CT hy = (CT)hyd;
  const CT hz = (CT)hzd;
  const CT rcut2 = (CT)(rcut*rcut);
  const CT inv_2sigmasq = (CT)(1.0/(2.0*sigma*sigma));
  const CT pref = (CT)pow(2.0*pi_dbl*sigma*sigma,-3.0/2.0);

  // Clear charge grid
  rho.clear();

  printf("CpuGaussCharge::spreadChargeToGrid, rcut=%lf nx=%d support=%d\n",
	 rcut,nx,calcSupportSize(rcut, boxx, boxy, boxz));

  /*
  // Spread charge on grid
  CT* rhodata = rho.getDataPointer();
  for (int i=0;i < numCoord;i++) {
    CT x = xyzq[i].x;
    CT y = xyzq[i].y;
    CT z = xyzq[i].z;
    CT q = xyzq[i].q*pref;
    int ix = (int)round(x*inv_hx);
    int iy = (int)round(y*inv_hy);
    int iz = (int)round(z*inv_hz);
    for (int tz=iz-nz;tz <= iz+nz;tz++) {
      CT zp = ((CT)tz)*hz;
      CT dz = z-zp;
      int pz = (tz + sizeZ) % sizeZ;
      for (int ty=iy-ny;ty <= iy+ny;ty++) {
	CT yp = ((CT)ty)*hy;
	CT dy = y-yp;
	int py = (ty + sizeY) % sizeY;
	for (int tx=ix-nx;tx <= ix+nx;tx++) {
	  CT xp = ((CT)tx)*hx;
	  CT dx = x-xp;
	  int px = (tx + sizeX) % sizeX;
	  CT r2 = dx*dx + dy*dy + dz*dz;
	  if (r2 < rcut2) {
	    int p = px + (py + pz*sizeY)*sizeX;
	    rhodata[p] += q*exp(-r2*inv_2sigmasq);
	  }
	}
      }
    }
  }
  */
  
  CT* wx = new CT[2*nx+1];
  CT* wy = new CT[2*ny+1];
  CT* wz = new CT[2*nz+1];
  
  // Spread charge on grid
  CT* rhodata = rho.getDataPointer();
  for (int i=0;i < numCoord;i++) {
    CT x = xyzq[i].x;
    CT y = xyzq[i].y;
    CT z = xyzq[i].z;
    CT q = xyzq[i].q*pref;
    int ixc = (int)round(x*inv_hx) - nx;
    int iyc = (int)round(y*inv_hy) - ny;
    int izc = (int)round(z*inv_hz) - nz;
    const CT dx = ixc*hx - x;
    const CT dy = iyc*hy - y;
    const CT dz = izc*hz - z;
    // Calculate 1d weights
    for (int j=-0;j <= 2*nx;j++) {
      CT x = dx + j*hx;
      wx[j] = expf(-x*x*inv_2sigmasq);
    }
    for (int j=-0;j <= 2*ny;j++) {
      CT y = dy + j*hy;
      wy[j] = expf(-y*y*inv_2sigmasq);
    }
    for (int j=-0;j <= 2*nz;j++) {
      CT z = dz + j*hz;
      wz[j] = q*expf(-z*z*inv_2sigmasq);
    }
    // Make sure (ixc, iyc, izc) are non-negative
    ixc = (ixc + sizeX) % sizeX;
    iyc = (iyc + sizeY) % sizeY;
    izc = (izc + sizeZ) % sizeZ;
    //
    for (int tz=0;tz <= 2*nz;tz++) {
      for (int ty=0;ty <= 2*ny;ty++) {
	for (int tx=0;tx <= 2*nx;tx++) {
	  int ix = (ixc + tx) % sizeX;
	  int iy = (iyc + ty) % sizeY;
	  int iz = (izc + tz) % sizeZ;
	  const int p = ix + (iy + iz*sizeY)*sizeX;
	  rhodata[p] += wx[tx]*wy[ty]*wz[tz];
	}
      }
    }
  }
  
  delete [] wx;
  delete [] wy;
  delete [] wz;
}

//
// On-grid charge spreading
//
template <typename AT, typename CT>
void CpuGaussCharge<AT, CT>::spreadChargeOnGrid(const double sigma, const double rcut,
						const double boxx, const double boxy, const double boxz,
						const CpuGrid<CT>& rho, CpuGrid<CT>& rhoS) {
  // Sanity checks
  assert(sizeX == rho.getSizeX());
  assert(sizeY == rho.getSizeY());
  assert(sizeZ == rho.getSizeZ());
  assert(sizeX == rhoS.getSizeX());
  assert(sizeY == rhoS.getSizeY());
  assert(sizeZ == rhoS.getSizeZ());
  
  const double hxd = boxx/(double)sizeX;
  const double hyd = boxy/(double)sizeY;
  const double hzd = boxz/(double)sizeZ;
  // NOTE: assumes lambdaSigma2 = lambdaSigma1
  //const double rSigma2d = this->lambdaSigma1*sqrt(2.0)*this->sigma2;
  const int nx = (int)ceil(rcut/hxd);
  const int ny = (int)ceil(rcut/hyd);
  const int nz = (int)ceil(rcut/hzd);
  const CT inv_2sigmasq = (CT)(1.0/(2.0*sigma*sigma));
  const CT pref = (CT)(pow(2.0*pi_dbl*sigma*sigma,-3.0/2.0)*hxd*hyd*hzd);
  const CT hxsq = (CT)(hxd*hxd);
  const CT hysq = (CT)(hyd*hyd);
  const CT hzsq = (CT)(hzd*hzd);

  if (rhoX == NULL) rhoX = new CpuGrid<CT>(sizeX, sizeY, sizeZ);
  if (rhoXY == NULL) rhoXY = new CpuGrid<CT>(sizeX, sizeY, sizeZ);

  printf("CpuGaussCharge::spreadChargeOnGrid, nx=%d\n",nx);
  
  // Pre-compute (Gx, Gy, Gz)
  // NOTE: prefix "pref" is included in Gz[]
  CT *Gx = new CT[2*nx+1];
  CT *Gy = new CT[2*ny+1];
  CT *Gz = new CT[2*nz+1];
  for (int ii=-nx;ii <= nx;ii++) {
    Gx[ii+nx] = exp(-hxsq*(CT)(ii*ii)*inv_2sigmasq);
  }
  for (int ii=-ny;ii <= ny;ii++) {
    Gy[ii+ny] = exp(-hysq*(CT)(ii*ii)*inv_2sigmasq);
  }
  for (int ii=-nz;ii <= nz;ii++) {
    Gz[ii+nz] = pref*exp(-hzsq*(CT)(ii*ii)*inv_2sigmasq);
  }
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();

  const CT* rhodata = rho.getDataPointer();
  CT* rhoXdata = rhoX->getDataPointer();
  CT* rhoXYdata = rhoXY->getDataPointer();
  CT* rhoSdata = rhoS.getDataPointer();

  // Clear charge grids
  rhoX->clear();
  rhoXY->clear();
  rhoS.clear();

  // Convolute
  for (int pz=0;pz < sizeZ;pz++) {
    for (int py=0;py < sizeY;py++) {
      int pyz = (py + pz*sizeY)*sizeX;
      for (int ix=0;ix < sizeX;ix++) {
	CT sum = (CT)0.0;
	for (int ii=-nx;ii <= nx;ii++) {
	  int px = (ii+ix+sizeX) % sizeX;
	  int p = px + pyz;
	  sum += Gx[ii+nx]*rhodata[p];
	}
	rhoXdata[ix + pyz] = sum;
	/*
	CT val = rhodata[ix + pyz];
	for (int ii=-nx;ii <= nx;ii++) {
	  int px = (ii+ix+sizeX) % sizeX;
	  int p = px + pyz;
	  rhoXdata[p] += Gx[ii+nx]*val;
	}
	*/
      }
    }
  }

  for (int pz=0;pz < sizeZ;pz++) {
    for (int px=0;px < sizeX;px++) {
      for (int iy=0;iy < sizeY;iy++) {
	CT sum = (CT)0.0;
	for (int ii=-ny;ii <= ny;ii++) {
	  int py = (ii+iy+sizeY) % sizeY;
	  int p = px + (py + pz*sizeY)*sizeX;
	  sum += Gy[ii+ny]*rhoXdata[p];
	}
	rhoXYdata[px + (iy + pz*sizeY)*sizeX] = sum;
	/*
	CT val = rhoXdata[px + (iy + pz*sizeY)*sizeX];
	for (int ii=-ny;ii <= ny;ii++) {
	  int py = (ii+iy+sizeY) % sizeY;
	  int p = px + (py + pz*sizeY)*sizeX;
	  rhoXYdata[p] += Gy[ii+ny]*val;
	}
	*/
      }
    }
  }

  for (int px=0;px < sizeX;px++) {
    for (int py=0;py < sizeY;py++) {
      for (int iz=0;iz < sizeZ;iz++) {
	CT sum = (CT)0.0;
	for (int ii=-nz;ii <= nz;ii++) {
	  int pz = (ii+iz+sizeZ) % sizeZ;
	  int p = px + (py + pz*sizeY)*sizeX;
	  sum += Gz[ii+nz]*rhoXYdata[p];
	}
	rhoSdata[px + (py + iz*sizeY)*sizeX] = sum;
	/*
	CT val = rhoXYdata[px + (py + iz*sizeY)*sizeX];
	for (int ii=-nz;ii <= nz;ii++) {
	  int pz = (ii+iz+sizeZ) % sizeZ;
	  int p = px + (py + pz*sizeY)*sizeX;
	  rhoSdata[p] += Gz[ii+nz]*val;
	}
	*/
      }
    }
  }

  delete [] Gx;
  delete [] Gy;
  delete [] Gz;
}

//
// Explicit instances of CpuGaussCharge
//
template class CpuGaussCharge<float, float>;
template class CpuGaussCharge<double, float>;
template class CpuGaussCharge<double, double>;
