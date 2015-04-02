#include <iostream>
#include <cassert>
#include "CpuFFTSolver.h"

//
// Class creator
//
template <typename T>
CpuFFTSolver<T>::CpuFFTSolver(CpuGreensFunc<T>& greensFunc, const int sizeX, const int sizeY, const int sizeZ) :
  greensFunc(greensFunc), sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ), rhoK((sizeX/2+1)*2, sizeY, sizeZ),
  plan_r2c_3d_created(false), plan_r2c_3d_in_pointer(NULL),
  plan_c2r_3d_created(false), plan_c2r_3d_out_pointer(NULL){
}

//
// Class destructor
//
template <typename T>
CpuFFTSolver<T>::~CpuFFTSolver() {
  if (plan_r2c_3d_created) fftw_destroy_plan(plan_r2c_3d);
  if (plan_c2r_3d_created) fftw_destroy_plan(plan_c2r_3d);
}

//
// Setup plans
//
template <typename T>
void CpuFFTSolver<T>::setupPlans(CpuGrid<T>& phi, const CpuGrid<T>& rho) {
  // R2C
  T* rhoPointer = const_cast<T *>(rho.getDataPointer());
  if (rhoPointer != plan_r2c_3d_in_pointer) {
    // Pointer has changed, re-create plan
    if (plan_r2c_3d_created) {
      fftw_destroy_plan(plan_r2c_3d);
      plan_r2c_3d_created = false;
    }
    if (sizeof(T) == 4) {
      std::cerr << "CpuFFTSolver::setupPlans, float type not yet implemented" << std::endl;
      exit(1);
    } else if (sizeof(T) == 8) {
      plan_r2c_3d = fftw_plan_dft_r2c_3d(sizeZ, sizeY, sizeX,
					 (double *)rhoPointer,
					 (fftw_complex *)rhoK.getDataPointer(),
					 FFTW_ESTIMATE);
    }
    plan_r2c_3d_created = true;
    plan_r2c_3d_in_pointer = rhoPointer;
  }
  // C2R
  T* phiPointer = phi.getDataPointer();
  if (phiPointer != plan_c2r_3d_out_pointer) {
    // Pointer has changed, re-create plan
    if (plan_c2r_3d_created) {
      fftw_destroy_plan(plan_c2r_3d);
      plan_c2r_3d_created = false;
    }
    if (sizeof(T) == 4) {
      std::cerr << "CpuFFTSolver::setupPlans, float type not yet implemented" << std::endl;
      exit(1);
    } else if (sizeof(T) == 8) {
      plan_c2r_3d = fftw_plan_dft_c2r_3d(sizeZ, sizeY, sizeX,
					 (fftw_complex *)rhoK.getDataPointer(),
					 (double *)phiPointer,
					 FFTW_ESTIMATE);
    }
    plan_c2r_3d_created = true;
    plan_c2r_3d_out_pointer = phiPointer;
  }  
}

//
// Solve!
//
template <typename T>
void CpuFFTSolver<T>::run(CpuGrid<T>& phi, const CpuGrid<T>& rho, const double sigma,
			  const double boxx, const double boxy, const double boxz) {
  // Sanity checks
  assert(sizeX == phi.getSizeX());
  assert(sizeY == phi.getSizeY());
  assert(sizeZ == phi.getSizeZ());
  assert(sizeX == rho.getSizeX());
  assert(sizeY == rho.getSizeY());
  assert(sizeZ == rho.getSizeZ());
  // Setup FFTW plan
  setupPlans(phi, rho);
  // Transform from Real to Complex: rho -> rhoK
  //rho.save("rho.txt");
  fftw_execute(plan_r2c_3d);
  //rhoK.save("rhoK1.txt");
  // Multiply rhoK by the Greens function
  greensFunc.mul(sizeX, sizeY, sizeZ, boxx, boxy, boxz, sigma, rhoK);
  //rhoK.save("rhoK2.txt");
  // Transform from Complex to Real: rhoK -> phi
  fftw_execute(plan_c2r_3d);
  //phi.save("phi.txt");
}

//
// Explicit instances of CpuFFTSolver
//
template class CpuFFTSolver<float>;
template class CpuFFTSolver<double>;
