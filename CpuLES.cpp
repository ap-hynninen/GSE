#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include "TypesGSE.h"
#include "CpuLES.h"

#define E_CURL_TYPE 2

//
// Class creator
//
template <typename AT, typename CT>
CpuLES<AT,CT>::CpuLES(const int sizeX, const int sizeY, const int sizeZ, const double sigma,
		      const double boxx, const double boxy, const double boxz,
		      const double ccelec) : 
  sigma(sigma), boxx(boxx), boxy(boxy), boxz(boxz), ccelec(ccelec),
  Ex(sizeX, sizeY, sizeZ), Ey(sizeX, sizeY, sizeZ), Ez(sizeX, sizeY, sizeZ),
  Bx(sizeX, sizeY, sizeZ), By(sizeX, sizeY, sizeZ), Bz(sizeX, sizeY, sizeZ),
  rhoTmp(sizeX, sizeY, sizeZ), rho(sizeX, sizeY, sizeZ), gaussCharge(sizeX, sizeY, sizeZ) {

  JxLen = 0;
  Jx = NULL;
  JyLen = 0;
  Jy = NULL;
  JzLen = 0;
  Jz = NULL;
}

//
// Class destructor
//
template <typename AT, typename CT>
CpuLES<AT,CT>::~CpuLES() {
  if (Jx != NULL) delete [] Jx;
  if (Jy != NULL) delete [] Jy;
  if (Jz != NULL) delete [] Jz;
}

//
// Calculate total energy
//
template <typename AT, typename CT>
double CpuLES<AT,CT>::calcTotalEnergy() {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const double hd = boxx/(double)sizeX;
  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();
  double Utot = 0.0;
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int n  = ix + (iy + iz*sizeY)*sizeX;
	int mi = (ix-1+sizeX)%sizeX + (iy + iz*sizeY)*sizeX;
	int pi = (ix+1)%sizeX + (iy + iz*sizeY)*sizeX;
	int mj = ix + ((iy-1+sizeY)%sizeY + iz*sizeY)*sizeX;
	int pj = ix + ((iy+1)%sizeY + iz*sizeY)*sizeX;
	int mk = ix + (iy + ((iz-1+sizeZ)%sizeZ)*sizeY)*sizeX;
	int pk = ix + (iy + ((iz+1)%sizeZ)*sizeY)*sizeX;

	double Ex_mi = ExData[mi];
	double Ex_nn = ExData[n];
	double Ex_pi = ExData[pi];

	double Ey_mj = EyData[mj];
	double Ey_nn = EyData[n];
	double Ey_pj = EyData[pj];

	double Ez_mk = EzData[mk];
	double Ez_nn = EzData[n];
	double Ez_pk = EzData[pk];

	Utot +=
	  (5.0/6.0)*(Ex_nn*Ex_nn) + (1.0/12.0)*Ex_nn*(Ex_mi + Ex_pi) +
	  (5.0/6.0)*(Ey_nn*Ey_nn) + (1.0/12.0)*Ey_nn*(Ey_mj + Ey_pj) +
	  (5.0/6.0)*(Ez_nn*Ez_nn) + (1.0/12.0)*Ez_nn*(Ez_mk + Ez_pk);
      }
    }
  }
  Utot *= ccelec*hd*hd*hd/(8.0*pi_dbl);
  return Utot;
}

//
// Calculates dipole sum = sum(qi*ri)
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::calcDipoleSum(const int numCoord, const xyzq_t<CT> *xyzq,
				  AT& dipSumX, AT& dipSumY, AT& dipSumZ) {
  dipSumX = (AT)0.0;
  dipSumY = (AT)0.0;
  dipSumZ = (AT)0.0;
  for (int i=0;i < numCoord;i++) {
    CT x = xyzq[i].x;
    CT y = xyzq[i].y;
    CT z = xyzq[i].z;
    CT q = xyzq[i].q;
    dipSumX += (AT)(q*x);
    dipSumY += (AT)(q*y);
    dipSumZ += (AT)(q*z);
  }
}

//
// Spread (interpolate) charge onto grid
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::spreadCharge1(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq) {
  gaussCharge.spreadChargeToGrid(sigma1, lambdaSigma1*sqrt(2.0)*sigma1, numCoord, xyzq,
				 boxx, boxy, boxz, rho);
}

//
// Spread charge on grid
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::spreadCharge2(const double sigma2, const double lambdaSigma1) {
  gaussCharge.spreadChargeOnGrid(sigma2, lambdaSigma1*sqrt(2.0)*sigma2,
				 boxx, boxy, boxz, rhoTmp, rho);
}

//
// Clear magnetic field
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::clearMagneticField() {
  Bx.clear();
  By.clear();
  Bz.clear();
}

//
// Sets electric field
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::setElectricField(const CpuGrid<CT>& ExIn, const CpuGrid<CT>& EyIn, const CpuGrid<CT>& EzIn) {
  Ex.copy(ExIn);
  Ey.copy(EyIn);
  Ez.copy(EzIn);
}

//
// Checks the validity of Gauss' Law. Returns the maximum error
//
template <typename AT, typename CT>
double CpuLES<AT,CT>::checkGaussLaw() {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const double hd = boxx/(double)sizeX;
  const CT fourpi_h = (CT)(4.0*pi_dbl*hd);
  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();
  CT* rhoData = rho.getDataPointer();
  double max_err = 0.0;
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int p = ix + (iy + iz*sizeY)*sizeX;
	int px = (ix-1+sizeX)%sizeX + (iy + iz*sizeY)*sizeX;
	int py = ix + ((iy-1+sizeY)%sizeY + iz*sizeY)*sizeX;
	int pz = ix + (iy + (iz-1+sizeZ)%sizeZ*sizeY)*sizeX;
	CT Ep = (ExData[p] - ExData[px]) + (EyData[p] - EyData[py]) + (EzData[p] - EzData[pz]);
	double err = fabs(rhoData[p]*fourpi_h - Ep);
	max_err = (err > max_err) ? err : max_err;
      }
    }
  }
  return max_err;
}

//
// Initialize electric field, from Joerg Rottler code (maxwell.cpp)
// NOT WORKING CORRECTLY YET
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::initElectricFieldJR() {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const int size = rho.getSize();

  const double hd = boxx/(double)sizeX;
  const CT fourpi_h = (CT)(4.0*pi_dbl*hd);

  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();
  const CT* rhoData = rho.getDataPointer();
  
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      CT dq = 0;
      for (int ix=0;ix < sizeX-1;ix++) {
	int p = ix + (iy + iz*sizeY)*sizeX;
	ExData[p] += rhoData[p]*fourpi_h + dq;
	dq = ExData[p];
      }
    }
  }

  for (int iz=0;iz < sizeZ;iz++) {
    CT dq = 0;
    for (int iy=0;iy < sizeY-1;iy++) {
      int p = (iy + iz*sizeY)*sizeX;
      EyData[p + sizeX-1] += rhoData[p + sizeX-1]*fourpi_h + ExData[p + sizeX-2] + dq;
      dq = EyData[p + sizeX-1];
    }
  }

  CT dq = 0;
  for (int iz=0;iz < sizeZ;iz++) {
    int p = iz*sizeY*sizeX;
    EzData[p + (sizeY-1)*sizeX + sizeX-1] += rhoData[p + (sizeY-1)*sizeX + sizeX-1]*fourpi_h +
      EyData[p + (sizeY-2)*sizeX + sizeX-1] +
      ExData[p + (sizeY-1)*sizeX + sizeX-2] + dq;
    dq = EzData[p + (sizeY-1)*sizeX + sizeX-1];
  }
  
}

//
// Initialize electric field
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::initElectricField() {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const int size = rho.getSize();
  //const double inv_hd = (double)sizeX/boxx;
  // NOTE: factor 4*pi is included in inv_h2
  //const CT inv_h2 = (CT)(4.0*pi_dbl*inv_hd*inv_hd);

  const double hd = boxx/(double)sizeX;
  const CT fourpi_h = (CT)(4.0*pi_dbl*hd);

  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();
  CT* rhoData = rho.getDataPointer();
  bool *ExSet = new bool[size];
  bool *EySet = new bool[size];
  bool *EzSet = new bool[size];
  for (int i=0;i < size;i++) {
    ExSet[i] = false;
    EySet[i] = false;
    EzSet[i] = false;
    ExData[i] = (CT)0.0;
    EyData[i] = (CT)0.0;
    EzData[i] = (CT)0.0;
  }
  int dix=1;
  int diy=1;
  int ix=0;
  int iy=0;
  int prev_ix=-1;
  int prev_iy=-1;
  int prev_iz=-1;
  int nvisit=0;
  for (int iz=0;iz < sizeZ;iz++) {
    for (;iy >=0 && iy < sizeY;iy+=diy) {
      for (;ix >= 0 && ix < sizeX;ix+=dix) {
	// ijk
	int n = ix + (iy + iz*sizeY)*sizeX;
	if (prev_ix != -1) {
	  int p = prev_ix + (prev_iy + prev_iz*sizeY)*sizeX;
	  int px = (prev_ix-1+sizeX)%sizeX + (prev_iy + prev_iz*sizeY)*sizeX;
	  int py = prev_ix + ((prev_iy-1+sizeY)%sizeY + prev_iz*sizeY)*sizeX;
	  int pz = prev_ix + (prev_iy + (prev_iz-1+sizeZ)%sizeZ*sizeY)*sizeX;
	  // Sanity check
	  if (ExSet[p] && EySet[p] && EzSet[p] && ExSet[px] && EySet[py] && EzSet[pz]) {
	    std::cout << "CpuLES::initElectricField, error: n=" << n
		      << " ix=" << ix << " iy=" << iy << " iz=" << iz << std::endl;
	    exit(1);
	  }
	  // Moved in +x direction
	  if (ix-prev_ix == 1) {
	    if (ExSet[p]) {
	      std::cout << "ExSet: n=" << n << std::endl
			<< "prev_ix,iy,iz = " << prev_ix << " " << prev_iy << " " << prev_iz << std::endl
			<< "ix,iy,iz      = " << ix << " " << iy << " " << iz << std::endl;
	      exit(1);
	    }
	    double Ep = EyData[p] + EzData[p] - ExData[px] - EyData[py] - EzData[pz];
	    ExData[p] = rhoData[p]*fourpi_h - Ep;
	  }
	  // Moved in -x direction
	  if (ix-prev_ix == -1) {
	    if (ExSet[px]) {
	      std::cout << "ExSet: n=" << n << std::endl
			<< "prev_ix,iy,iz = " << prev_ix << " " << prev_iy << " " << prev_iz << std::endl
			<< "ix,iy,iz      = " << ix << " " << iy << " " << iz << std::endl;
	      exit(1);
	    }
	    double Ep = ExData[p] + EyData[p] + EzData[p] - EyData[py] - EzData[pz];
	    ExData[px] = -(rhoData[p]*fourpi_h - Ep);
	  }
	  // Moved in +y direction
	  if (iy-prev_iy == 1) {
	    if (EySet[p]) {
	      std::cout << "EySet: n=" << n << std::endl
			<< "prev_ix,iy,iz = " << prev_ix << " " << prev_iy << " " << prev_iz << std::endl
			<< "ix,iy,iz      = " << ix << " " << iy << " " << iz << std::endl;
	      exit(1);
	    }
	    double Ep = ExData[p] + EzData[p] - ExData[px] - EyData[py] - EzData[pz];
	    EyData[p] = rhoData[p]*fourpi_h - Ep;
	  }
	  // Moved in -y direction
	  if (iy-prev_iy == -1) {
	    if (EySet[py]) {
	      std::cout << "EySet: n=" << n << std::endl
			<< "prev_ix,iy,iz = " << prev_ix << " " << prev_iy << " " << prev_iz << std::endl
			<< "ix,iy,iz      = " << ix << " " << iy << " " << iz << std::endl;
	      exit(1);
	    }
	    double Ep = ExData[p] + EyData[p] + EzData[p] - ExData[px] - EzData[pz];
	    EyData[py] = -(rhoData[p]*fourpi_h - Ep);
	  }
	  // Moved in +z direction
	  if (iz-prev_iz == 1) {
	    if (EzSet[p]) {
	      std::cout << "EzSet: n=" << n << std::endl
			<< "prev_ix,iy,iz = " << prev_ix << " " << prev_iy << " " << prev_iz << std::endl
			<< "ix,iy,iz      = " << ix << " " << iy << " " << iz << std::endl;
	      exit(1);
	    }
	    double Ep = ExData[p] + EyData[p] - ExData[px] - EyData[py] - EzData[pz];
	    EzData[p] = rhoData[p]*fourpi_h - Ep;
	  }
	  // Moved in -z direction
	  if (iz-prev_iz == -1) {
	    std::cout << "Moves in -z direction not possible" << std::endl;
	    exit(1);
	  }

	  double Ep = ExData[p] + EyData[p] + EzData[p] - ExData[px] - EyData[py] - EzData[pz];
	  double err = rhoData[p]*fourpi_h - Ep;
	  if (fabs(err) > 1.0e-10) {
	    printf("err=%lf\n",err);
	    exit(1);
	  }
	  
	  ExSet[p] = true;
	  EySet[p] = true;
	  EzSet[p] = true;
	  ExSet[px] = true;
	  EySet[py] = true;
	  EzSet[pz] = true;
	}
	prev_ix = ix;
	prev_iy = iy;
	prev_iz = iz;
	nvisit++;
      }
      // Flip x direction
      dix = -dix;
      ix += dix;
    }
    // Flip y direction
    diy = -diy;
    iy += diy;
  }
  //std::cout << "CpuLES::initElectricField, nvisit = " << nvisit << std::endl;
  int n_not_set=0;
  for (int iz=0;iz < sizeZ;iz++) {
    for (iy=0;iy < sizeY;iy++) {
      for (ix=0;ix < sizeX;ix++) {
	int n = ix + (iy + iz*sizeY)*sizeX;
	if (!ExSet[n]) {
	  std::cout << "Ex not set: " << " ix=" << ix << " iy=" << iy << " iz=" << iz << std::endl;
	  n_not_set++;
	}
	if (!EySet[n]) {
	  std::cout << "Ey not set: " << " ix=" << ix << " iy=" << iy << " iz=" << iz << std::endl;
	  n_not_set++;
	}
	if (!EzSet[n]) {
	  std::cout << "Ez not set: " << " ix=" << ix << " iy=" << iy << " iz=" << iz << std::endl;
	  n_not_set++;
	}
      }
    }
  }

  delete [] ExSet;
  delete [] EySet;
  delete [] EzSet;
}

//
// Interpolate force
//
template <typename AT, typename CT>
void CpuLES<AT, CT>::interpolateForce(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq,
				      AT *forceX, AT *forceY, AT *forceZ) {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const double rcut = lambdaSigma1*sqrt(2.0)*sigma1;
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
  const CT inv_2sigmasq = (CT)(1.0/(2.0*sigma1*sigma1));

  const CT* ExData = Ex.getDataPointer();
  const CT* EyData = Ey.getDataPointer();
  const CT* EzData = Ez.getDataPointer();

  printf("CpuLES::interpolateForceJ, rcut=%lf nx,ny,nz=%d %d %d\n",rcut,nx,ny,nz);

  for (int i=0;i < numCoord;i++) {
    CT x = xyzq[i].x;
    CT y = xyzq[i].y;
    CT z = xyzq[i].z;
    int ix = (int)round(x*inv_hx);
    int iy = (int)round(y*inv_hy);
    int iz = (int)round(z*inv_hz);
    AT fx = (AT)0.0;
    AT fy = (AT)0.0;
    AT fz = (AT)0.0;
    int j0 = i*(2*nx+1)*(2*ny+1)*(2*nz+1);
    for (int tz=-nz;tz <= nz;tz++) {
      int oz = (iz+tz + sizeZ) % sizeZ;
      for (int ty=-ny;ty <= ny;ty++) {
	int oy = (iy+ty + sizeY) % sizeY;
	for (int tx=-nx;tx <= nx;tx++) {
	  int ox = (ix+tx + sizeX) % sizeX;

	  int n = ox + (oy + oz*sizeY)*sizeX;
	  int mx = (ox-1+sizeX)%sizeX + (oy + oz*sizeY)*sizeX;
	  int px = (ox+1)%sizeX + (oy + oz*sizeY)*sizeX;
	  int my = ox + ((oy-1+sizeY)%sizeY + oz*sizeY)*sizeX;
	  int py = ox + ((oy+1)%sizeY + oz*sizeY)*sizeX;
	  int mz = ox + (oy + ((oz-1+sizeZ)%sizeZ)*sizeY)*sizeX;
	  int pz = ox + (oy + ((oz+1)%sizeZ)*sizeY)*sizeX;

	  int jn = j0 + (tx+nx) + ((ty+ny) + (tz+nz)*(2*ny+1))*(2*nx+1);
	  int jmx = j0 + (tx-1+nx) + ((ty+ny) + (tz+nz)*(2*ny+1))*(2*nx+1);
	  int jpx = j0 + (tx+1+nx) + ((ty+ny) + (tz+nz)*(2*ny+1))*(2*nx+1);

	  int jmy = j0 + (tx+nx) + ((ty-1+ny) + (tz+nz)*(2*ny+1))*(2*nx+1);
	  int jpy = j0 + (tx+nx) + ((ty+1+ny) + (tz+nz)*(2*ny+1))*(2*nx+1);

	  int jmz = j0 + (tx+nx) + ((ty+ny) + (tz-1+nz)*(2*ny+1))*(2*nx+1);
	  int jpz = j0 + (tx+nx) + ((ty+ny) + (tz+1+nz)*(2*ny+1))*(2*nx+1);

	  CT Exm = ExData[mx];
	  CT Exn = ExData[n];
	  CT Exp = ExData[px];

	  CT Eym = EyData[my];
	  CT Eyn = EyData[n];
	  CT Eyp = EyData[py];

	  CT Ezm = EzData[mz];
	  CT Ezn = EzData[n];
	  CT Ezp = EzData[pz];

	  CT Jxm = (tx-1 >= -nx) ? Jx[jmx] : 0;
	  CT Jxn = Jx[jn];
	  CT Jxp = (tx+1 <=  nx) ? Jx[jpx] : 0;

	  CT Jym = (ty-1 >= -ny) ? Jy[jmy] : 0;
	  CT Jyn = Jy[jn];
	  CT Jyp = (ty+1 <=  ny) ? Jy[jpy] : 0;

	  CT Jzm = (tz-1 >= -nz) ? Jz[jmz] : 0;
	  CT Jzn = Jz[jn];
	  CT Jzp = (tz+1 <=  nz) ? Jz[jpz] : 0;

	  fx += Jxn*(5.0/6.0*Exn + 1.0/24.0*(Exm + Exp)) + 1.0/24.0*(Jxm + Jxp)*Exn;
	  fy += Jyn*(5.0/6.0*Eyn + 1.0/24.0*(Eym + Eyp)) + 1.0/24.0*(Jym + Jyp)*Eyn;
	  fz += Jzn*(5.0/6.0*Ezn + 1.0/24.0*(Ezm + Ezp)) + 1.0/24.0*(Jzm + Jzp)*Ezn;
	}
      }
    }
    forceX[i] = -fx*ccelec*hx*hx*hx;
    forceY[i] = -fy*ccelec*hy*hy*hy;
    forceZ[i] = -fz*ccelec*hz*hz*hz;
  }
 
  // Dipole part
  AT dipSumX, dipSumY, dipSumZ;
  gaussCharge.calcDipoleSum(boxx, boxy, boxz, rho, dipSumX, dipSumY, dipSumZ);
  AT EdipX = 4.0*pi_dbl/(boxx*boxy*boxz)*dipSumX*ccelec;
  AT EdipY = 4.0*pi_dbl/(boxx*boxy*boxz)*dipSumY*ccelec;
  AT EdipZ = 4.0*pi_dbl/(boxx*boxy*boxz)*dipSumZ*ccelec;
  printf("Edip: %e %e %e\n",EdipX,EdipY,EdipZ);
  for (int i=0;i < numCoord;i++) {
    CT q = xyzq[i].q;
    forceX[i] += q*EdipX;
    forceY[i] += q*EdipY;
    forceZ[i] += q*EdipZ;
  }
}

//
// Implements charge fluctuation. This updates electric fields and calculates current
// that is later used to compute forces
// xyzq_p = previous coordinates
// xyzq_c = current coordinates
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::chargeFluctuation(const double sigma, const double lambdaSigma,
				      const int numCoord, const xyzq_t<CT> *xyzq_p, const xyzq_t<CT> *xyzq_c) {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const double hxd = boxx/(double)sizeX;
  const double hyd = boxy/(double)sizeY;
  const double hzd = boxz/(double)sizeZ;
  const CT inv_hx = (CT)((double)sizeX/boxx);
  const CT inv_hy = (CT)((double)sizeY/boxy);
  const CT inv_hz = (CT)((double)sizeZ/boxz);
  const CT hx = (CT)(boxx/(double)sizeX);
  const CT hy = (CT)(boxy/(double)sizeY);
  const CT hz = (CT)(boxz/(double)sizeZ);
  const double rcut = lambdaSigma*sqrt(2.0)*sigma;
  const int nx = (int)ceil(rcut/hxd);
  const int ny = (int)ceil(rcut/hyd);
  const int nz = (int)ceil(rcut/hzd);
  const CT sigmasq = sigma*sigma;
  const CT pref = (CT)pow(2.0*pi_dbl*sigmasq,-1.0/2.0);
  assert(hx == hy);
  assert(hx == hz);
  const CT inv_sigma2sq = (CT)(hx*hx/(2.0*sigmasq));
  const CT inv_sigmasq = (CT)(hx*hx/sigmasq);

  int JreqLen = numCoord*(2*nx+1)*(2*ny+1)*(2*nz+1);
  if (Jx != NULL) {
    if (JxLen < JreqLen) {
      delete [] Jx;
      Jx = NULL;
    }
  }
  if (Jx == NULL) {
    JxLen = JreqLen;
    Jx = new CT[JxLen];
  }
  
  if (Jy != NULL) {
    if (JyLen < JreqLen) {
      delete [] Jy;
      Jy = NULL;
    }
  }
  if (Jy == NULL) {
    JyLen = JreqLen;
    Jy = new CT[JyLen];
  }

  if (Jz != NULL) {
    if (JzLen < JreqLen) {
      delete [] Jz;
      Jz = NULL;
    }
  }
  if (Jz == NULL) {
    JzLen = JreqLen;
    Jz = new CT[JzLen];
  }

  CT* wpx = new CT[2*nx+1];
  CT* wpy = new CT[2*ny+1];
  CT* wpz = new CT[2*nz+1];
  CT* dwpx = new CT[2*nx+1];
  CT* dwpy = new CT[2*ny+1];
  CT* dwpz = new CT[2*nz+1];

  CT* wmx = new CT[2*nx+1];
  CT* wmy = new CT[2*ny+1];
  CT* wmz = new CT[2*nz+1];

  CT* wcx = new CT[2*nx+1];
  CT* wcy = new CT[2*ny+1];
  CT* wcz = new CT[2*nz+1];

  CT* delq1x = new CT[2*nx+1];
  CT* delq1y = new CT[2*ny+1];
  CT* delq1z = new CT[2*nz+1];
  CT* jlinkx = new CT[2*nx+1];
  CT* jlinky = new CT[2*ny+1];
  CT* jlinkz = new CT[2*nz+1];

  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();

  printf("CpuLES::chargeFluctuation, rcut=%lf nx,ny,nz=%d %d %d\n",rcut,nx,ny,nz);
  
  for (int i=0;i < numCoord;i++) {
    // Current position
    CT xc = xyzq_c[i].x;
    CT yc = xyzq_c[i].y;
    CT zc = xyzq_c[i].z;
    int ixc = (int)round(xc*inv_hx);
    int iyc = (int)round(yc*inv_hy);
    int izc = (int)round(zc*inv_hz);

    // Previous position
    CT xp = xyzq_p[i].x;
    CT yp = xyzq_p[i].y;
    CT zp = xyzq_p[i].z;
    int ixp = (int)round(xp*inv_hx);
    int iyp = (int)round(yp*inv_hy);
    int izp = (int)round(zp*inv_hz);

    // Midpoint
    CT xm = 0.5*(xc+xp);
    CT ym = 0.5*(xc+xp);
    CT zm = 0.5*(xc+xp);
    int ixm = (int)round(xm*inv_hx);
    int iym = (int)round(ym*inv_hy);
    int izm = (int)round(zm*inv_hz);

    CT dx, dy, dz;
    // Previous position
    dx = (CT)ixp - xp*inv_hx;
    dy = (CT)iyp - yp*inv_hy;
    dz = (CT)izp - zp*inv_hz;
    for (int j=-nx;j <= nx;j++) wpx[j+nx] = pref*exp(-(dx+j)*(dx+j)*inv_sigma2sq);
    for (int j=-ny;j <= ny;j++) wpy[j+ny] = pref*exp(-(dy+j)*(dy+j)*inv_sigma2sq);
    for (int j=-nz;j <= nz;j++) wpz[j+nz] = pref*exp(-(dz+j)*(dz+j)*inv_sigma2sq);
    for (int j=-nx;j <= nx;j++) dwpx[j+nx] = wpx[j+nx]*((-dx-j)*inv_sigmasq);
    for (int j=-ny;j <= ny;j++) dwpy[j+ny] = wpy[j+ny]*((-dy-j)*inv_sigmasq);
    for (int j=-nz;j <= nz;j++) dwpz[j+nz] = wpz[j+nz]*((-dz-j)*inv_sigmasq);
   
    // Mid position
    dx = (CT)ixm - xm*inv_hx;
    dy = (CT)iym - ym*inv_hy;
    dz = (CT)izm - zm*inv_hz;
    for (int j=-nx;j <= nx;j++) wmx[j+nx] = pref*exp(-(dx+j)*(dx+j)*inv_sigma2sq);
    for (int j=-ny;j <= ny;j++) wmy[j+ny] = pref*exp(-(dy+j)*(dy+j)*inv_sigma2sq);
    for (int j=-nz;j <= nz;j++) wmz[j+nz] = pref*exp(-(dz+j)*(dz+j)*inv_sigma2sq);

    // Current position
    dx = (CT)ixc - xc*inv_hx;
    dy = (CT)iyc - yc*inv_hy;
    dz = (CT)izc - zc*inv_hz;
    for (int j=-nx;j <= nx;j++) wcx[j+nx] = pref*exp(-(dx+j)*(dx+j)*inv_sigma2sq);
    for (int j=-ny;j <= ny;j++) wcy[j+ny] = pref*exp(-(dy+j)*(dy+j)*inv_sigma2sq);
    for (int j=-nz;j <= nz;j++) wcz[j+nz] = pref*exp(-(dz+j)*(dz+j)*inv_sigma2sq);

    delq1x[0] = wcx[0] - wpx[0];
    jlinkx[0] = dwpx[0];
    for (int j=-nx+1;j <= nx;j++) {
      delq1x[j+nx] = delq1x[j-1+nx] + (wcx[j+nx] - wpx[j+nx]);
      jlinkx[j+nx] = jlinkx[j-1+nx] + dwpx[j+nx];
    }

    delq1y[0] = wcy[0] - wpy[0];
    jlinky[0] = dwpy[0];
    for (int j=-ny+1;j <= ny;j++) {
      delq1y[j+ny] = delq1y[j-1+ny] + (wcy[j+ny] - wpy[j+ny]);
      jlinky[j+ny] = jlinky[j-1+ny] + dwpy[j+ny];
    }

    delq1z[0] = wcz[0] - wpz[0];
    jlinkz[0] = dwpz[0];
    for (int j=-nz+1;j <= nz;j++) {
      delq1z[j+nz] = delq1z[j-1+nz] + (wcz[j+nz] - wpz[j+nz]);
      jlinkz[j+nz] = jlinkz[j-1+nz] + dwpz[j+nz];
    }

    int j0 = i*(2*nx+1)*(2*ny+1)*(2*nz+1);
    for (int tz=-nz;tz<=nz;tz++) {
      int iz = izc + tz;
      for (int ty=-ny;ty<=ny;ty++) {
	int iy = iyc + ty;
	for (int tx=-nx;tx<=nx;tx++) {
	  int ix = ixc + tx;
	  int p = ((ix+sizeX)%sizeX) + ( ((iy+sizeY)%sizeY) + ((iz+sizeZ)%sizeZ)*sizeY )*sizeX;
	  CT tmp1 = xyzq_c[i].q*wpy[ty+ny]*wpz[tz+nz];
	  CT tmp2 = xyzq_c[i].q*wcx[tx+nx]*wpz[tz+nz];
	  CT tmp3 = xyzq_c[i].q*wcx[tx+nx]*wcy[ty+ny];
	  ExData[p] += tmp1*delq1x[tx+nx];
	  EyData[p] += tmp2*delq1y[ty+ny];
	  EzData[p] += tmp3*delq1z[tz+nz];
	  int jp = j0 + (tx+nx) + ((ty+ny) + (tz+nz)*(2*ny+1))*(2*nx+1);
	  Jx[jp] = -tmp1*jlinkx[tx+nx];
	  Jy[jp] = -tmp2*jlinky[ty+ny];
	  Jz[jp] = -tmp3*jlinkz[tz+nz];
	}
      }
    }
    
  }

  delete [] wpx;
  delete [] wpy;
  delete [] wpz;
  delete [] dwpx;
  delete [] dwpy;
  delete [] dwpz;

  delete [] wmx;
  delete [] wmy;
  delete [] wmz;

  delete [] wcx;
  delete [] wcy;
  delete [] wcz;

  delete [] delq1x;
  delete [] delq1y;
  delete [] delq1z;
  delete [] jlinkx;
  delete [] jlinky;
  delete [] jlinkz;
}

//
// Integrate Maxwell equations
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::integrate(const CT c, const CT dt, const CT gamma2) {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const int size = rho.getSize();
  const CT inv_h = 1.0;//(CT)((double)sizeX/boxx);
  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();
  CT* BxData = Bx.getDataPointer();
  CT* ByData = By.getDataPointer();
  CT* BzData = Bz.getDataPointer();
  CT one_min_gamma2_dt = (CT)1.0 - gamma2*dt;
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	// ijk
	int n  = ix + (iy + iz*sizeY)*sizeX;
	CT curlx, curly, curlz;
	curl(E_CURL_TYPE, ix, iy, iz, sizeX, sizeY, sizeZ, inv_h, ExData, EyData, EzData, curlx, curly, curlz);
	BxData[n] = one_min_gamma2_dt*BxData[n] - dt*c*c*curlx;
	ByData[n] = one_min_gamma2_dt*ByData[n] - dt*c*c*curly;
	BzData[n] = one_min_gamma2_dt*BzData[n] - dt*c*c*curlz;
      }
    }
  }

  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int n  = ix + (iy + iz*sizeY)*sizeX;
	CT curlx, curly, curlz;
	curl(-1, ix, iy, iz, sizeX, sizeY, sizeZ, inv_h, BxData, ByData, BzData, curlx, curly, curlz);
	ExData[n] += dt*curlx;
	EyData[n] += dt*curly;
	EzData[n] += dt*curlz;
      }
    }
  }

}

//
// Returns curl of D
// curlType = -1: (D[i]   - D[i-1])/h
// curlType =  1: (D[i+1] - D[i])/h
// curlType =  2: ( (5/6)*(D[i+1] - D[i]) + (1/24)*(D[i] - D[i-1]) + (1/24)*(D[i+2] - D[i+1]))
//
template <typename AT, typename CT>
CT CpuLES<AT,CT>::curl(const int curlType, const int ix, const int iy, const int iz,
		       const int sizeX, const int sizeY, const int sizeZ,
		       const CT inv_h,
		       const CT* DxData, const CT* DyData, const CT* DzData,
		       CT& valx, CT& valy, CT& valz) {
  if (curlType == -1) {
    // ijk (n=neutral position)
    int n  = ix + (iy + iz*sizeY)*sizeX;
    // i-1jk (m=minus 1 position)
    int mi = (ix-1+sizeX)%sizeX + (iy + iz*sizeY)*sizeX;
    // ij-1k
    int mj = ix + ((iy-1+sizeY)%sizeY + iz*sizeY)*sizeX;
    // ijk-1
    int mk = ix + (iy + ((iz-1+sizeZ)%sizeZ)*sizeY)*sizeX;
    valx = inv_h*DzData[n] - DzData[mj] - (DyData[n] - DyData[mk]);
    valy = inv_h*DxData[n] - DxData[mk] - (DzData[n] - DzData[mi]);
    valz = inv_h*DyData[n] - DyData[mi] - (DxData[n] - DxData[mj]);
  } else if (curlType == 1) {
    // ijk (n=neutral position)
    int n  = ix + (iy + iz*sizeY)*sizeX;
    // i+1jk (m=minus 1 position)
    int pi = (ix+1)%sizeX + (iy + iz*sizeY)*sizeX;
    // ij+1k
    int pj = ix + ((iy+1)%sizeY + iz*sizeY)*sizeX;
    // ijk+1
    int pk = ix + (iy + ((iz+1)%sizeZ)*sizeY)*sizeX;
    valx = inv_h*DzData[pj] - DzData[n] - (DyData[pk] - DyData[n]);
    valy = inv_h*DxData[pk] - DxData[n] - (DzData[pi] - DzData[n]);
    valz = inv_h*DyData[pi] - DyData[n] - (DxData[pj] - DxData[n]);
  } else if (curlType == 2) {
    // ijk (n=neutral position)
    int n  = ix + (iy + iz*sizeY)*sizeX;
    // i-1jk
    int mi = (ix-1+sizeX)%sizeX + (iy + iz*sizeY)*sizeX;
    // i+1jk
    int pi = (ix+1)%sizeX + (iy + iz*sizeY)*sizeX;
    // i+2jk
    int pi2 = (ix+2)%sizeX + (iy + iz*sizeY)*sizeX;
    // ij-1k
    int mj = ix + ((iy-1+sizeY)%sizeY + iz*sizeY)*sizeX;
    // ij+1k
    int pj = ix + ((iy+1)%sizeY + iz*sizeY)*sizeX;
    // ij+2k
    int pj2 = ix + ((iy+2)%sizeY + iz*sizeY)*sizeX;
    // ijk-1
    int mk = ix + (iy + ((iz-1+sizeZ)%sizeZ)*sizeY)*sizeX;
    // ijk+1
    int pk = ix + (iy + ((iz+1)%sizeZ)*sizeY)*sizeX;
    // ijk+2
    int pk2 = ix + (iy + ((iz+2)%sizeZ)*sizeY)*sizeX;
    valx = inv_h*((CT)(5.0/6.0)*(DzData[pj]  - DzData[n] -  (DyData[pk]  - DyData[n])) +
		  (CT)(1.0/24.0)*(DzData[n]   - DzData[mj] - (DyData[n]   - DyData[mk])) + 
		  (CT)(1.0/24.0)*(DzData[pj2] - DzData[pj] - (DyData[pk2] - DyData[pk])));
    valy = inv_h*((CT)(5.0/6.0)*(DxData[pk]  - DxData[n] -  (DzData[pi]  - DzData[n])) +
		  (CT)(1.0/24.0)*(DxData[n]   - DxData[mk] - (DzData[n]   - DzData[mi])) + 
		  (CT)(1.0/24.0)*(DxData[pk2] - DxData[pk] - (DzData[pi2] - DzData[pi])));
    valz = inv_h*((CT)(5.0/6.0)*(DyData[pi]  - DyData[n] -  (DxData[pj]  - DxData[n])) +
		  (CT)(1.0/24.0)*(DyData[n]   - DyData[mi] - (DxData[n]   - DxData[mj])) + 
		  (CT)(1.0/24.0)*(DyData[pi2] - DyData[pi] - (DxData[pj2] - DxData[pj])));
  }
}

//
// Returns maximum norm of curl
//
template <typename AT, typename CT>
CT CpuLES<AT,CT>::maxCurl(const int curlType, const CpuGrid<CT>& Dx, const CpuGrid<CT>& Dy, const CpuGrid<CT>& Dz) {
  const int sizeX = Dx.getSizeX();
  const int sizeY = Dx.getSizeY();
  const int sizeZ = Dx.getSizeZ();
  const CT* DxData = Dx.getDataPointer();
  const CT* DyData = Dy.getDataPointer();
  const CT* DzData = Dz.getDataPointer();
  const CT inv_h = (CT)((double)sizeX/boxx);
  CT maxvalx = (CT)0.0;
  CT maxvaly = (CT)0.0;
  CT maxvalz = (CT)0.0;
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	CT valx, valy, valz;
	curl(curlType, ix, iy, iz, sizeX, sizeY, sizeZ, inv_h, DxData, DyData, DzData, valx, valy, valz);
	valx = fabs(valx);
	valy = fabs(valy);
	valz = fabs(valz);
	maxvalx = (valx > maxvalx) ? valx : maxvalx;
	maxvaly = (valy > maxvaly) ? valy : maxvaly;
	maxvalz = (valz > maxvalz) ? valz : maxvalz;
      }
    }
  }
  //printf("maxval = %e %e %e\n",maxvalx,maxvaly,maxvalz);
  CT maxval = (maxvalx > maxvaly) ? maxvalx : maxvaly;
  maxval = (maxvalz > maxval) ? maxvalz : maxval;
  return maxval;
}

//
// Returns maximum norm of curl(B)
//
template <typename AT, typename CT>
CT CpuLES<AT,CT>::maxCurlB() {
  return maxCurl(-1, Bx, By, Bz);
}

//
// Returns maximum norm of curl(E)
//
template <typename AT, typename CT>
CT CpuLES<AT,CT>::maxCurlE() {
  return maxCurl(E_CURL_TYPE, Ex, Ey, Ez);
}

//
// Calculates dipole term energy
// NOTE: erf(L/(sqrt(2)*S)) term is dropped since it's usually close to 1.0
//
template <typename AT, typename CT>
double CpuLES<AT,CT>::calcDipoleEnergy(const int numCoord, const xyzq_t<CT>* xyzq) {
  AT dipSumX, dipSumY, dipSumZ;
  gaussCharge.calcDipoleSum(boxx, boxy, boxz, rho, dipSumX, dipSumY, dipSumZ);
  printf("dipSum = %lf %lf %lf\n",dipSumX, dipSumY, dipSumZ);
  return 2.0*pi_dbl/(boxx*boxy*boxz)*(dipSumX*dipSumX + dipSumY*dipSumY + dipSumZ*dipSumZ)*ccelec;
}

//
// Explicit instances of this class
//
template class CpuLES<double,float>;
template class CpuLES<double,double>;
