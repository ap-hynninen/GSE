#include <cassert>
#include <cmath>
#include <iostream>
#include "CpuMultiGridSolver.h"

//
// Class creator
//
template <typename T>
CpuMultiGridSolver<T>::CpuMultiGridSolver(const int sizeFinestX, const int sizeFinestY, const int sizeFinestZ) {
  assert(sizeFinestX == sizeFinestY);
  assert(sizeFinestX == sizeFinestZ);
  
  // Goal is to coarse down to size 4
  numLevel = 1;//(int)log2((double)sizeFinestX/16.0) + 1;
  rhoVec.reserve(numLevel);
  phiVec.reserve(numLevel);
  resVec.reserve(numLevel);

  // level = 0:          finest
  // level = numLevel-1  coarsest
  int sizeX = sizeFinestX;
  int sizeY = sizeFinestY;
  int sizeZ = sizeFinestZ;
  rhoVec.push_back(NULL);
  phiVec.push_back(NULL);
  resVec.push_back(new CpuGrid<T>(sizeX, sizeY, sizeZ));
  for (int indLevel=1;indLevel < numLevel;indLevel++) {
    sizeX /= 2;
    sizeY /= 2;
    sizeZ /= 2;
    rhoVec.push_back(new CpuGrid<T>(sizeX, sizeY, sizeZ));
    phiVec.push_back(new CpuGrid<T>(sizeX, sizeY, sizeZ));
    resVec.push_back(new CpuGrid<T>(sizeX, sizeY, sizeZ));
  }

  numPreRelax = 10;
  numPostRelax = 10;  
}

//
// Class destructor
//
template <typename T>
CpuMultiGridSolver<T>::~CpuMultiGridSolver() {
  for (int i=0;i < numLevel;i++) {
    if (rhoVec.at(i) != NULL) delete rhoVec.at(i);
    if (phiVec.at(i) != NULL) delete phiVec.at(i);
    delete resVec.at(i);
  }
}			  

//
// Run solver for Poisson equation: D^2 phi = rho
// rho = charge distribution i.e. source term
// phi = electrostatic potential i.e. result
//
template <typename T>
void CpuMultiGridSolver<T>::run(CpuGrid<T>& phi, const CpuGrid<T>& rho,
				const double boxx, const double boxy, const double boxz) {

  for (int icycle=0;icycle < 1;icycle++) {
  
    // Go down V: fine -> coarse
    for (int i=0;i < numLevel-1;i++) {
      const CpuGrid<T>* rhoFine = (i == 0) ? &rho : rhoVec.at(i);
      CpuGrid<T>* phiFine = (i == 0) ? &phi : phiVec.at(i);
      // Relax (pre-smoothing)
      //int laplaceRep = (i == 0) ? HERMITIAN : SEVEN_POINT;
      int laplaceRep = SEVEN_POINT;
      for (int j=0;j < numPreRelax;j++) relaxRBGS(*phiFine, *rhoFine, boxx, boxy, boxz, laplaceRep);
      // Calculate residual
      CpuGrid<T>* resFine = resVec.at(i);
      calculateResidual(*phiFine, *rhoFine, *resFine, boxx, boxy, boxz, laplaceRep);
      printf("down V residual: %e\n",resFine->maxAbsValue());
      // Move to next coarser level
      CpuGrid<T>* rhoCoarse = rhoVec.at(i+1);
      CpuGrid<T>* phiCoarse = phiVec.at(i+1);
      fineToCoarse(*resFine, *rhoCoarse);
      // Clear next initial guess
      phiCoarse->clear();
    }
  
    // Find solution at the coarsest grid
    const CpuGrid<T>* rhoCoarsest = (numLevel <= 1) ? &rho : rhoVec.at(numLevel-1);
    CpuGrid<T>* phiCoarsest = (numLevel <= 1) ? &phi : phiVec.at(numLevel-1);

    int laplaceRepFinest = HERMITIAN;
    for (int j=0;j < 200;j++) {
      relaxRBGS(*phiCoarsest, *rhoCoarsest, boxx, boxy, boxz, laplaceRepFinest);
    }
    phiCoarsest->save("phiCoarsest.txt");
    rhoCoarsest->save("rhoCoarsest.txt");

    calculateResidual(*phiCoarsest, *rhoCoarsest, *resVec.at(numLevel-1), boxx, boxy, boxz, laplaceRepFinest);
    resVec.at(numLevel-1)->save("resCoarsest.txt");
    printf("residual: %e\n",resVec.at(numLevel-1)->maxAbsValue());
  
    // Go up V: coarse -> fine
    for (int i=numLevel-1;i >= 1;i--) {
      // Move to finer level, store temporary result into resFine
      const CpuGrid<T>* phiCoarse = phiVec.at(i);
      CpuGrid<T>* tmp = resVec.at(i-1);
      coarseToFine(*phiCoarse, *tmp);
      // Add: phiFine += P[phiCoarse]
      CpuGrid<T>* phiFine = (i == 1) ? &phi : phiVec.at(i-1); 
      phiFine->add(*tmp);
      // Relax (post-smoothing)
      const CpuGrid<T>* rhoFine = (i == 1) ? &rho : rhoVec.at(i-1); 
      //int laplaceRep = (i == 1) ? HERMITIAN : SEVEN_POINT;
      int laplaceRep = SEVEN_POINT;
      for (int j=0;j < numPostRelax;j++) relaxRBGS(*phiFine, *rhoFine, boxx, boxy, boxz, laplaceRep);
      //
      calculateResidual(*phiFine, *rhoFine, *tmp, boxx, boxy, boxz, laplaceRep);
      printf("up V residual: %e\n",tmp->maxAbsValue());
    }
  }
}

//
// Calculate residual d = L*phi - rho
//
template <typename T>
void CpuMultiGridSolver<T>::calculateResidual(const CpuGrid<T>& phi, const CpuGrid<T>& rho, CpuGrid<T>& res,
					      const double boxx, const double boxy, const double boxz,
					      const int laplaceRep) {
  // Sanity checks
  assert(phi.getSizeX() == rho.getSizeX());
  assert(phi.getSizeY() == rho.getSizeY());
  assert(phi.getSizeZ() == rho.getSizeZ());
  assert(phi.getSizeX() == res.getSizeX());
  assert(phi.getSizeY() == res.getSizeY());
  assert(phi.getSizeZ() == res.getSizeZ());
  
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  T* resData = res.getDataPointer();

  double inv_hx2d = (double)sizeX/boxx;
  double inv_hy2d = (double)sizeY/boxy;
  double inv_hz2d = (double)sizeZ/boxz;
  inv_hx2d *= inv_hx2d;
  inv_hy2d *= inv_hy2d;
  inv_hz2d *= inv_hz2d;
  T inv_hx2 = (T)inv_hx2d;
  T inv_hy2 = (T)inv_hy2d;
  T inv_hz2 = (T)inv_hz2d;

  T a_A;
  T inv_a_A;
  T a_B;
  T b_B;
  T bx_A, by_A, bz_A;
  T cxy_A, cxz_A, cyz_A;
  if (laplaceRep == HERMITIAN) {
    double sum_inv = inv_hx2d + inv_hy2d + inv_hz2d;
    double a_Ad = 4.0/3.0*sum_inv;
    a_A = (T)a_Ad;
    inv_a_A = (T)(1.0/a_Ad);
    a_B = (T)(-0.5);
    b_B = (T)(-1.0/12.0);
    bx_A = (T)(-5.0/6.0*inv_hx2d + sum_inv/6.0);
    by_A = (T)(-5.0/6.0*inv_hy2d + sum_inv/6.0);
    bz_A = (T)(-5.0/6.0*inv_hz2d + sum_inv/6.0);
    cxy_A = (T)(-1.0/12.0*(inv_hx2d + inv_hy2d));
    cxz_A = (T)(-1.0/12.0*(inv_hx2d + inv_hz2d));
    cyz_A = (T)(-1.0/12.0*(inv_hy2d + inv_hz2d));
  }

  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int p = ix + (iy + iz*sizeY)*sizeX;
	if (laplaceRep == SEVEN_POINT) {
	  T twoPhiVal = ((T)2.0)*phi.getDataValue(ix, iy, iz);
	  resData[p] = ((phi.getDataValue(ix-1, iy, iz) - twoPhiVal + phi.getDataValue(ix+1, iy, iz))*inv_hx2 +
			(phi.getDataValue(ix, iy-1, iz) - twoPhiVal + phi.getDataValue(ix, iy+1, iz))*inv_hy2 +
			(phi.getDataValue(ix, iy, iz-1) - twoPhiVal + phi.getDataValue(ix, iy, iz+1))*inv_hz2 -
			rho.getDataValue(ix, iy, iz));
	} else {
	  T B_rho = a_B*rho.getDataValue(ix, iy, iz) +
	    b_B*(rho.getDataValue(ix-1, iy, iz) + rho.getDataValue(ix+1, iy, iz)) +
	    b_B*(rho.getDataValue(ix, iy-1, iz) + rho.getDataValue(ix, iy+1, iz)) +
	    b_B*(rho.getDataValue(ix, iy, iz-1) + rho.getDataValue(ix, iy, iz+1));
	  T bA_phi = bx_A*(phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz)) +
	    by_A*(phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz)) +
	    bz_A*(phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1));
	  T cA_phi = cxy_A*(phi.getDataValue(ix-1, iy-1, iz) + phi.getDataValue(ix-1, iy+1, iz) +
			    phi.getDataValue(ix+1, iy-1, iz) + phi.getDataValue(ix+1, iy+1, iz)) +
	    cxz_A*(phi.getDataValue(ix-1, iy, iz-1) + phi.getDataValue(ix-1, iy, iz+1) +
		   phi.getDataValue(ix+1, iy, iz-1) + phi.getDataValue(ix+1, iy, iz+1)) +
	    cyz_A*(phi.getDataValue(ix, iy-1, iz-1) + phi.getDataValue(ix, iy-1, iz+1) +
		   phi.getDataValue(ix, iy+1, iz-1) + phi.getDataValue(ix, iy+1, iz+1));
	  resData[p] = a_A*phi.getDataValue(ix, iy, iz) + bA_phi + cA_phi -B_rho;
	}
      }
    }
  }
}

//
// Relax grid using red-black Gauss-Seidel method
//
template <typename T>
void CpuMultiGridSolver<T>::relaxRBGS(CpuGrid<T>& phi, const CpuGrid<T>& rho,
				      const double boxx, const double boxy, const double boxz,
				      const int laplaceRep) {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  T* phiData = phi.getDataPointer();

  // Sanity checks
  assert(rho.getSizeX() == phi.getSizeX());
  assert(rho.getSizeY() == phi.getSizeY());
  assert(rho.getSizeZ() == phi.getSizeZ());
  assert(sizeX % 2 == 0);
  assert(sizeY % 2 == 0);
  assert(sizeZ % 2 == 0);

  double inv_hx2d = (double)sizeX/boxx;
  double inv_hy2d = (double)sizeY/boxy;
  double inv_hz2d = (double)sizeZ/boxz;
  inv_hx2d *= inv_hx2d;
  inv_hy2d *= inv_hy2d;
  inv_hz2d *= inv_hz2d;
  T fac = (T)(1.0/(2.0*inv_hx2d + 2.0*inv_hy2d + 2.0*inv_hz2d));
  T inv_hx2 = (T)inv_hx2d;
  T inv_hy2 = (T)inv_hy2d;
  T inv_hz2 = (T)inv_hz2d;

  T inv_a_A;
  T a_B;
  T b_B;
  T bx_A, by_A, bz_A;
  T cxy_A, cxz_A, cyz_A;
  if (laplaceRep == HERMITIAN) {
    double sum_inv = inv_hx2d + inv_hy2d + inv_hz2d;
    double a_A = 4.0/3.0*sum_inv;
    inv_a_A = (T)(1.0/a_A);
    a_B = (T)(-0.5);
    b_B = (T)(-1.0/12.0);
    bx_A = (T)(-5.0/6.0*inv_hx2d + sum_inv/6.0);
    by_A = (T)(-5.0/6.0*inv_hy2d + sum_inv/6.0);
    bz_A = (T)(-5.0/6.0*inv_hz2d + sum_inv/6.0);
    cxy_A = (T)(-1.0/12.0*(inv_hx2d + inv_hy2d));
    cxz_A = (T)(-1.0/12.0*(inv_hx2d + inv_hz2d));
    cyz_A = (T)(-1.0/12.0*(inv_hy2d + inv_hz2d));
  }

  int iz0 = 0;
  for (int ipass=0;ipass < 2;ipass++) {
    int iy0 = iz0;
    for (int iz=0;iz < sizeZ;iz++) {
      int ix0 = iy0;
      for (int iy=0;iy < sizeY;iy++) {
	for (int ix=ix0;ix < sizeX;ix+=2) {
	  int p = ix + (iy + iz*sizeY)*sizeX;
	  if (laplaceRep == SEVEN_POINT) {
	    phiData[p] = ((phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz))*inv_hx2 +
			  (phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz))*inv_hy2 +
			  (phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1))*inv_hz2 -
			  rho.getDataValue(ix, iy, iz))*fac;
	  } else {
	    T B_rho = a_B*rho.getDataValue(ix, iy, iz) +
	      b_B*(rho.getDataValue(ix-1, iy, iz) + rho.getDataValue(ix+1, iy, iz)) +
	      b_B*(rho.getDataValue(ix, iy-1, iz) + rho.getDataValue(ix, iy+1, iz)) +
	      b_B*(rho.getDataValue(ix, iy, iz-1) + rho.getDataValue(ix, iy, iz+1));
	    T bA_phi = bx_A*(phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz)) +
	      by_A*(phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz)) +
	      bz_A*(phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1));
	    T cA_phi = cxy_A*(phi.getDataValue(ix-1, iy-1, iz) + phi.getDataValue(ix-1, iy+1, iz) +
			      phi.getDataValue(ix+1, iy-1, iz) + phi.getDataValue(ix+1, iy+1, iz)) +
	      cxz_A*(phi.getDataValue(ix-1, iy, iz-1) + phi.getDataValue(ix-1, iy, iz+1) +
		     phi.getDataValue(ix+1, iy, iz-1) + phi.getDataValue(ix+1, iy, iz+1)) +
	      cyz_A*(phi.getDataValue(ix, iy-1, iz-1) + phi.getDataValue(ix, iy-1, iz+1) +
		     phi.getDataValue(ix, iy+1, iz-1) + phi.getDataValue(ix, iy+1, iz+1));
	    phiData[p] = inv_a_A*(B_rho - bA_phi - cA_phi);
	  }
	}
	ix0 = (ix0 + 1) % 2;
      }
      iy0 = (iy0 + 1) % 2;
    }
    iz0 = (iz0 + 1) % 2;
  }

  /*
  for (int ipass=0;ipass < 2;ipass++) {
    for (int iz=0;iz < sizeZ;iz++) {
      for (int iy=0;iy < sizeY;iy++) {
	for (int ix=0;ix < sizeX;ix++) {
	  if ((ix + iy + iz) % 2 == ipass) {
	    int p = ix + (iy + iz*sizeY)*sizeX;
	    phiData[p] = ((phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz))*inv_hx2 +
			  (phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz))*inv_hy2 +
			  (phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1))*inv_hz2 -
			  rho.getDataValue(ix, iy, iz))*fac;
	  }
	}
      }
    }
  }
  */
  
}

//
// Relax grid using Gauss-Seidel method
//
template <typename T>
void CpuMultiGridSolver<T>::relaxGS(CpuGrid<T>& phi, const CpuGrid<T>& rho,
				    const double boxx, const double boxy, const double boxz,
				    const int laplaceRep) {
  const int sizeX = rho.getSizeX();
  const int sizeY = rho.getSizeY();
  const int sizeZ = rho.getSizeZ();
  T* phiData = phi.getDataPointer();

  // Sanity checks
  assert(rho.getSizeX() == phi.getSizeX());
  assert(rho.getSizeY() == phi.getSizeY());
  assert(rho.getSizeZ() == phi.getSizeZ());
  assert(sizeX % 2 == 0);
  assert(sizeY % 2 == 0);
  assert(sizeZ % 2 == 0);

  double inv_hx2d = (double)sizeX/boxx;
  double inv_hy2d = (double)sizeY/boxy;
  double inv_hz2d = (double)sizeZ/boxz;
  inv_hx2d *= inv_hx2d;
  inv_hy2d *= inv_hy2d;
  inv_hz2d *= inv_hz2d;
  T fac = (T)(1.0/(2.0*inv_hx2d + 2.0*inv_hy2d + 2.0*inv_hz2d));
  T inv_hx2 = (T)inv_hx2d;
  T inv_hy2 = (T)inv_hy2d;
  T inv_hz2 = (T)inv_hz2d;

  T inv_a_A;
  T a_B;
  T b_B;
  T bx_A, by_A, bz_A;
  T cxy_A, cxz_A, cyz_A;
  if (laplaceRep == HERMITIAN) {
    double sum_inv = inv_hx2d + inv_hy2d + inv_hz2d;
    double a_A = 4.0/3.0*sum_inv;
    inv_a_A = (T)(1.0/a_A);
    a_B = (T)(-0.5);
    b_B = (T)(-1.0/12.0);
    bx_A = (T)(-5.0/6.0*inv_hx2d + sum_inv/6.0);
    by_A = (T)(-5.0/6.0*inv_hy2d + sum_inv/6.0);
    bz_A = (T)(-5.0/6.0*inv_hz2d + sum_inv/6.0);
    cxy_A = (T)(-1.0/12.0*(inv_hx2d + inv_hy2d));
    cxz_A = (T)(-1.0/12.0*(inv_hx2d + inv_hz2d));
    cyz_A = (T)(-1.0/12.0*(inv_hy2d + inv_hz2d));
  }
  
  for (int iz=0;iz < sizeZ;iz++) {
    for (int iy=0;iy < sizeY;iy++) {
      for (int ix=0;ix < sizeX;ix++) {
	int p = ix + (iy + iz*sizeY)*sizeX;
	if (laplaceRep == SEVEN_POINT) {
	  phiData[p] = ((phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz))*inv_hx2 +
			(phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz))*inv_hy2 +
			(phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1))*inv_hz2 -
			rho.getDataValue(ix, iy, iz))*fac;
	} else {
	  T B_rho = a_B + b_B*(rho.getDataValue(ix-1, iy, iz) + rho.getDataValue(ix+1, iy, iz)) +
	    b_B*(rho.getDataValue(ix, iy-1, iz) + rho.getDataValue(ix, iy+1, iz)) + 
	    b_B*(rho.getDataValue(ix, iy, iz-1) + rho.getDataValue(ix, iy, iz+1));
	  T bA_phi = bx_A*(phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz)) +
	    by_A*(phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz)) +
	    bz_A*(phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1));
	  T cA_phi = cxy_A*(phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz) +
			    phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz)) +
	    cxz_A*(phi.getDataValue(ix-1, iy, iz) + phi.getDataValue(ix+1, iy, iz) +
		   phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1)) +
	    cyz_A*(phi.getDataValue(ix, iy-1, iz) + phi.getDataValue(ix, iy+1, iz) +
		   phi.getDataValue(ix, iy, iz-1) + phi.getDataValue(ix, iy, iz+1));
	  phiData[p] = inv_a_A*(B_rho - bA_phi - cA_phi);
	}
      }
    }
  }
}

//
// Returns the finite diffence Laplace operator of size as a matrix of size [size X size]
//
template <typename T>
void CpuMultiGridSolver<T>::buildLaplace(const CpuGrid<T>& phi, T* L, const double boxx, const double boxy, const double boxz) {
  const int sizeX = phi.getSizeX();
  const int sizeY = phi.getSizeY();
  const int sizeZ = phi.getSizeZ();
  const int size = sizeX*sizeY*sizeZ;
  double inv_hx2d = (double)sizeX/boxx;
  double inv_hy2d = (double)sizeY/boxy;
  double inv_hz2d = (double)sizeZ/boxz;
  inv_hx2d *= inv_hx2d;
  inv_hy2d *= inv_hy2d;
  inv_hz2d *= inv_hz2d;
  T inv_hx2 = (T)inv_hx2d;
  T inv_hy2 = (T)inv_hy2d;
  T inv_hz2 = (T)inv_hz2d;
  
  for (int i=0;i < size*size;i++) L[i] = (T)0;

  for (int z=0;z < sizeZ;z++) {
    for (int y=0;y < sizeY;y++) {
      for (int x=0;x < sizeX;x++) {
	int row = (x + (y + z*sizeY)*sizeX)*size;
	// double derivative w.r.t x
	L[phi.getPos(x-1, y, z) + row] += inv_hx2;
	L[phi.getPos(x,   y, z) + row] += (T)(-2.0)*inv_hx2;
	L[phi.getPos(x+1, y, z) + row] += inv_hx2;
	// double derivative w.r.t y
	L[phi.getPos(x, y-1, z) + row] += inv_hy2;
	L[phi.getPos(x, y,   z) + row] += (T)(-2.0)*inv_hy2;
	L[phi.getPos(x, y+1, z) + row] += inv_hy2;
	// double derivative w.r.t z
	L[phi.getPos(x, y, z-1) + row] += inv_hz2;
	L[phi.getPos(x, y, z  ) + row] += (T)(-2.0)*inv_hz2;
	L[phi.getPos(x, y, z+1) + row] += inv_hz2;
      }
    }
  }
  
}

//
// Exact solution at the coarsest grid
//
template <typename T>
void CpuMultiGridSolver<T>::solveSmall(CpuGrid<T>& phi, const CpuGrid<T>& rho) {
  // Sanity checks
  assert(phi.getSizeX() == rho.getSizeX());
  assert(phi.getSizeY() == rho.getSizeY());
  assert(phi.getSizeZ() == rho.getSizeZ());
  assert(phi.getSizeX() == 4);
  assert(phi.getSizeY() == 4);
  assert(phi.getSizeZ() == 4);

  
}

//
// Move from fine -> coarse grid
//
template <typename T>
void CpuMultiGridSolver<T>::fineToCoarse(const CpuGrid<T>& rhoFine, CpuGrid<T>& rhoCoarse) {
  const int sizeCoarseX = rhoCoarse.getSizeX();
  const int sizeCoarseY = rhoCoarse.getSizeY();
  const int sizeCoarseZ = rhoCoarse.getSizeZ();
  T* rhoCoarseData = rhoCoarse.getDataPointer();
  
  const int sizeFineX = rhoFine.getSizeX();
  const int sizeFineY = rhoFine.getSizeY();
  const int sizeFineZ = rhoFine.getSizeZ();
  const int sizeFine  = rhoFine.getSize();
  const T* rhoFineData = rhoFine.getDataPointer();

  // Sanity check
  assert(sizeFineX == 2*sizeCoarseX);
  assert(sizeFineY == 2*sizeCoarseY);
  assert(sizeFineZ == 2*sizeCoarseZ);

  fprintf(stderr,"fineToCoarse: %d %d %d -> %d %d %d\n",sizeFineX,sizeFineY,sizeFineZ,sizeCoarseX,sizeCoarseY,sizeCoarseZ);
  
  // Easy "Write once" algorithm
  for (int iz=0;iz < sizeCoarseZ;iz++) {
    for (int iy=0;iy < sizeCoarseY;iy++) {
      for (int ix=0;ix < sizeCoarseX;ix++) {
	T val = (T)0;
	for (int tz=-1;tz <= 1;tz++) {
	  T facZ = (tz == 0) ? (T)1 : (T)0.5;
	  int pz = (2*iz + tz + sizeFineZ) % sizeFineZ;
	  for (int ty=-1;ty <= 1;ty++) {
	    T facYZ = facZ*((ty == 0) ? (T)1 : (T)0.5);
	    int py = (2*iy + ty + sizeFineY) % sizeFineY;
	    for (int tx=-1;tx <= 1;tx++) {
	      T facXYZ = facYZ*((tx == 0) ? (T)1 : (T)0.5);
	      int px = (2*ix + tx + sizeFineX) % sizeFineX;
	      int p = px + (py + pz*sizeFineY)*sizeFineX;
	      val += facXYZ*rhoFineData[p];
	    }
	  }
	}
	val *= (T)0.25;
	int i = ix + (iy + iz*sizeCoarseY)*sizeCoarseX;
	rhoCoarseData[i] = val;
      }
    }
  }

}

//
// Move from coarse -> fine grid
//
template <typename T>
void CpuMultiGridSolver<T>::coarseToFine(const CpuGrid<T>& rhoCoarse, CpuGrid<T>& rhoFine) {
  const int sizeCoarseX = rhoCoarse.getSizeX();
  const int sizeCoarseY = rhoCoarse.getSizeY();
  const int sizeCoarseZ = rhoCoarse.getSizeZ();
  const T* rhoCoarseData = rhoCoarse.getDataPointer();
  
  const int sizeFineX = rhoFine.getSizeX();
  const int sizeFineY = rhoFine.getSizeY();
  const int sizeFineZ = rhoFine.getSizeZ();
  const int sizeFine  = rhoFine.getSize();
  T* rhoFineData = rhoFine.getDataPointer();

  // Sanity check
  assert(sizeFineX == 2*sizeCoarseX);
  assert(sizeFineY == 2*sizeCoarseY);
  assert(sizeFineZ == 2*sizeCoarseZ);

  fprintf(stderr,"coarseToFine: %d %d %d -> %d %d %d\n",sizeCoarseX,sizeCoarseY,sizeCoarseZ,sizeFineX,sizeFineY,sizeFineZ);

  /*
  // "write once" -algorithm
  for (int iz=0;iz < sizeFineZ;iz++) {
    for (int iy=0;iy < sizeFineY;iy++) {
      for (int ix=0;ix < sizeFineX;ix++) {
	rhoFineData[p] = 
      }
    }
  }
  */
  
  // Easy "Read once" -algorithm
  for (int i=0;i < sizeFine;i++) rhoFineData[i] = (T)0;
  for (int iz=0;iz < sizeCoarseZ;iz++) {
    for (int iy=0;iy < sizeCoarseY;iy++) {
      for (int ix=0;ix < sizeCoarseX;ix++) {
	T val = rhoCoarse.getDataValue(ix, iy, iz);
	for (int tz=-1;tz <= 1;tz++) {
	  T facZ = (tz == 0) ? (T)1 : (T)0.5;
	  int pz = (2*iz + tz + sizeFineZ) % sizeFineZ;
	  for (int ty=-1;ty <= 1;ty++) {
	    T facYZ = facZ*((ty == 0) ? (T)1 : (T)0.5);
	    int py = (2*iy + ty + sizeFineY) % sizeFineY;
	    for (int tx=-1;tx <= 1;tx++) {
	      T facXYZ = facYZ*((tx == 0) ? (T)1 : (T)0.5);
	      int px = (2*ix + tx + sizeFineX) % sizeFineX;
	      int p = px + (py + pz*sizeFineY)*sizeFineX;
	      rhoFineData[p] += facXYZ*val;
	    }
	  }
	}
      }
    }
  }
  
}

//
// Explicit instances of CpuMultiGridSolver
//
template class CpuMultiGridSolver<float>;
template class CpuMultiGridSolver<double>;
