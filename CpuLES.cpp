#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include "TypesGSE.h"
#include "CpuLES.h"

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
}

//
// Class destructor
//
template <typename AT, typename CT>
CpuLES<AT,CT>::~CpuLES() {}

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
  Utot *= ccelec*hd*hd*hd/(8.0*pi);
  return Utot;
}

//
// Interpolate force
//
template <typename AT, typename CT>
void CpuLES<AT, CT>::interpolateForceVW(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq,
					AT *forceX, AT *forceY, AT *forceZ) {

  gaussCharge.interpolateForce(sigma1, lambdaSigma1*sqrt(2.0)*sigma1, numCoord, xyzq,
			       boxx, boxy, boxz, Ex, Ey, Ez, forceX, forceY, forceZ);
  for (int i=0;i < numCoord;i++) {
    forceX[i] *= ccelec;
    forceY[i] *= ccelec;
    forceZ[i] *= ccelec;
  }
  AT dipSumX, dipSumY, dipSumZ;
  gaussCharge.calcDipoleSum(boxx, boxy, boxz, rho, dipSumX, dipSumY, dipSumZ);
  AT EdipX = 4.0*pi/(boxx*boxy*boxz)*dipSumX;
  AT EdipY = 4.0*pi/(boxx*boxy*boxz)*dipSumY;
  AT EdipZ = 4.0*pi/(boxx*boxy*boxz)*dipSumZ;
  for (int i=0;i < numCoord;i++) {
    CT q = (CT)(xyzq[i].q*ccelec);
    //printf("%lf\n",-q*EdipX);
    //forceX[i] -= q*EdipX;
    //forceY[i] -= q*EdipY;
    //forceZ[i] -= q*EdipZ;
  }

}

//
// Interpolate force
//
template <typename AT, typename CT>
void CpuLES<AT, CT>::interpolateForceEF(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq,
					AT *forceX, AT *forceY, AT *forceZ) {

  AT* ExPart = new AT[numCoord];
  AT* EyPart = new AT[numCoord];
  AT* EzPart = new AT[numCoord];

  interpolateElectricField(sigma1, lambdaSigma1, numCoord, xyzq, ExPart, EyPart, EzPart);
  for (int i=0;i < numCoord;i++) {
    CT q = xyzq[i].q*ccelec;
    forceX[i] = ExPart[i]*q;
    forceY[i] = EyPart[i]*q;
    forceZ[i] = EzPart[i]*q;
  }
  delete [] ExPart;
  delete [] EyPart;
  delete [] EzPart;

  // Dipole part
  AT dipSumX, dipSumY, dipSumZ;
  gaussCharge.calcDipoleSum(boxx, boxy, boxz, rho, dipSumX, dipSumY, dipSumZ);
  AT EdipX = 4.0*pi/(boxx*boxy*boxz)*dipSumX;
  AT EdipY = 4.0*pi/(boxx*boxy*boxz)*dipSumY;
  AT EdipZ = 4.0*pi/(boxx*boxy*boxz)*dipSumZ;
  
  for (int i=0;i < numCoord;i++) {
    CT q = xyzq[i].q*ccelec;
    //forceX[i] += q*EdipX;
    //forceY[i] += q*EdipY;
    //forceZ[i] += q*EdipZ;
  }

  printf("LES(dip): %e | %e\n",xyzq[0].q*ccelec*EdipX,xyzq[1].q*ccelec*EdipX);
}

//
// Interpolate force from mesh-based potential
//
template <typename AT, typename CT>
void CpuLES<AT, CT>::interpolateElectricField(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq,
					      AT *ExPart, AT *EyPart, AT *EzPart) {

  FILE* handle = fopen("E_LES.txt","wt");
  gaussCharge.interpolateElectricField(sigma1, lambdaSigma1*sqrt(2.0)*sigma1, numCoord, xyzq,
				       boxx, boxy, boxz, Ex, Ey, Ez, ExPart, EyPart, EzPart, handle);
  fclose(handle);

  /*
  // Calculate dipole term
  AT dipSumX, dipSumY, dipSumZ;
  calcDipoleSum(numCoord, xyzq, dipSumX, dipSumY, dipSumZ);
  //printf("dipSum = %e %e %e\n",dipSumX, dipSumY, dipSumZ);
  dipSumX *= (AT)(-2.0*sqrt(pi)/(boxx*boxy*boxz));
  dipSumY *= (AT)(-2.0*sqrt(pi)/(boxx*boxy*boxz));
  dipSumZ *= (AT)(-2.0*sqrt(pi)/(boxx*boxy*boxz));

  printf("dipSum = %e %e %e\n",dipSumX,dipSumY,dipSumZ);
  
  // Subtract dipole term
  for (int i=0;i < numCoord;i++) {
    ExPart[i] -= dipSumX;
    EyPart[i] -= dipSumY;
    EzPart[i] -= dipSumZ;
  }
  */
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
    /*
    while (x < (CT)0.0) x += (CT)boxx;
    while (x > (CT)boxx) x -= (CT)boxx;
    while (y < (CT)0.0) y += (CT)boxy;
    while (y > (CT)boxy) y -= (CT)boxy;
    while (z < (CT)0.0) z += (CT)boxz;
    while (z > (CT)boxz) z -= (CT)boxz;
    */
    CT q = xyzq[i].q;
    dipSumX += (AT)(q*x);
    dipSumY += (AT)(q*y);
    dipSumZ += (AT)(q*z);
  }
  /*
  while (dipSumX > (CT)boxx) dipSumX -= (CT)boxx;
  while (dipSumX < -(CT)boxx) dipSumX += (CT)boxx;
  while (dipSumY > (CT)boxy) dipSumY -= (CT)boxy;
  while (dipSumY < -(CT)boxy) dipSumY += (CT)boxy;
  while (dipSumZ > (CT)boxz) dipSumZ -= (CT)boxz;
  while (dipSumZ < -(CT)boxz) dipSumZ += (CT)boxz;
  */
}

//
// Spread (interpolate) charge onto grid
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::spreadCharge1(const double sigma1, const double lambdaSigma1, const int numCoord, const xyzq_t<CT> *xyzq) {
  gaussCharge.spreadChargeToGrid(sigma1, lambdaSigma1*sqrt(2.0)*sigma1, numCoord, xyzq,
				 boxx, boxy, boxz, rho);
  //const int sizeX = rho.getSizeX();
  //const CT h = (CT)(boxx/(double)sizeX);
  //printf("h*h*h = %lf\n",h*h*h);
  //rho.scale(h*h*h);
  //gaussCharge.spreadChargeToGrid(sigma, 3.0*sqrt(2.0)*sigma, numCoord, xyzq,
  //				 boxx, boxy, boxz, rho);
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
  const CT fourpi_h = (CT)(4.0*pi*hd);
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
// Initialize electric field
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::initElectricFieldHolm() {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  const int size = rho.getSize();
  const double inv_hd = (double)sizeX/boxx;
  // NOTE: factor 4*pi is included in inv_h2
  const CT inv_h2 = (CT)(4.0*pi*inv_hd*inv_hd);
  CT* ExData = Ex.getDataPointer();
  CT* EyData = Ey.getDataPointer();
  CT* EzData = Ez.getDataPointer();
  const CT* rhoData = rho.getDataPointer();
  bool *ExSet = new bool[size];
  bool *EySet = new bool[size];
  bool *EzSet = new bool[size];
  for (int i=0;i < size;i++) {
    ExData[i] = (CT)0.0;
    EyData[i] = (CT)0.0;
    EzData[i] = (CT)0.0;
  }
  for (int iz=0;iz < sizeZ;iz++) {
    // Plane
    CT Qz = (CT)0.0;
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int n = ix + (iy + iz*sizeY)*sizeX;
	Qz += rhoData[n];
      }
    }
    Qz /= (CT)(sizeX*sizeY);
    CT Qplane = Qz*inv_h2;
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int n = ix + (iy + iz*sizeY)*sizeX;
	int np = ix + (iy + ((iz+1)%sizeZ)*sizeY)*sizeX;
	EzData[np] = EzData[n] + Qplane;
      }
    }
    // Line
    for (int iy=0;iy < sizeY;iy++) {
      CT Qy = (CT)0.0;
      for (int ix=0;ix < sizeX;ix++) {
	int n = ix + (iy + iz*sizeY)*sizeX;
	Qy += rhoData[n];
      }
      Qy /= (CT)(sizeX);
      CT Qline = (Qy-Qz)*inv_h2;
      for (int ix=0;ix < sizeX;ix++) {
	int n = ix + (iy + iz*sizeY)*sizeX;
	int np = ix + ((iy+1)%sizeY + iz*sizeY)*sizeX;
	EyData[np] = EyData[n] + Qline;
      }
      // Vertex
      for (int ix=0;ix < sizeX;ix++) {
	int n = ix + (iy + iz*sizeY)*sizeX;
	int np = (ix+1)%sizeX + (iy + iz*sizeY)*sizeX;
	CT Qx = rhoData[n];
	ExData[np] = ExData[n] + (Qx-Qy)*inv_h2;
      }
    }
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
  //const CT inv_h2 = (CT)(4.0*pi*inv_hd*inv_hd);

  const double hd = boxx/(double)sizeX;
  const CT fourpi_h = (CT)(4.0*pi*hd);

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

  //std::cout << "CpuLES::initElectricField, n_not_set = " << n_not_set << std::endl;
 
  delete [] ExSet;
  delete [] EySet;
  delete [] EzSet;
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
	curl(2, ix, iy, iz, sizeX, sizeY, sizeZ, inv_h, ExData, EyData, EzData, curlx, curly, curlz);
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
// curlType =  2: ((D[i+1] - D[i])*10 + (D[i] - D[i-1]) + (D[i+2] - D[i+1]))/(12*h)
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
    valx = inv_h*( (CT)10.0*
		   (DzData[pj]  - DzData[n] -  (DyData[pk]  - DyData[n])) +
		   (DzData[n]   - DzData[mj] - (DyData[n]   - DyData[mk])) + 
		   (DzData[pj2] - DzData[pj] - (DyData[pk2] - DyData[pk]) ) )/(CT)(12.0);
    valy = inv_h*( (CT)10.0*
		   (DxData[pk]  - DxData[n] -  (DzData[pi]  - DzData[n])) +
		   (DxData[n]   - DxData[mk] - (DzData[n]   - DzData[mi])) + 
		   (DxData[pk2] - DxData[pk] - (DzData[pi2] - DzData[pi]) ) )/(CT)(12.0);
    valz = inv_h*( (CT)10.0*
		   (DyData[pi]  - DyData[n] -  (DxData[pj]  - DxData[n])) +
		   (DyData[n]   - DyData[mi] - (DxData[n]   - DxData[mj])) + 
		   (DyData[pi2] - DyData[pi] - (DxData[pj2] - DxData[pj]) ) )/(CT)(12.0);
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
  return maxCurl(2, Ex, Ey, Ez);
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
  return 2.0*pi/(boxx*boxy*boxz)*(dipSumX*dipSumX + dipSumY*dipSumY + dipSumZ*dipSumZ)*ccelec;
}

/*
//
// Save electric field in file
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::saveElectricField(const char* filename) {
  std::ofstream file(filename);
  if (file.is_open()) {
    for (int i=0;i < nx*sizeY*nz;i++) {
      file << Ex[i] << " " << Ey[i] << " " << Ez[i] << std::endl;
    }
  } else {
    std::cerr<<"Error opening file "<<filename<<std::endl;
    exit(1);
  }
}

//
// Save charge in file
//
template <typename AT, typename CT>
void CpuLES<AT,CT>::saveCharge(const char* filename) {
  std::ofstream file(filename);
  if (file.is_open()) {
    for (int i=0;i < nx*ny*nz;i++) {
      file << rho[i] << std::endl;
    }
  } else {
    std::cerr<<"Error opening file "<<filename<<std::endl;
    exit(1);
  }
}
*/

//
// Explicit instances of this class
//
template class CpuLES<double,float>;
template class CpuLES<double,double>;
