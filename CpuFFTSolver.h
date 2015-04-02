#ifndef CPUFFTSOLVER_H
#define CPUFFTSOLVER_H
//
// FFT solver for Poisson equation
// (c) Antti-Pekka Hynninen, aphynninen@hotmail.com
//
#include <fftw3.h>
#include "CpuGrid.h"
#include "CpuGreensFunc.h"

template <typename T>
class CpuFFTSolver {
private:
  //--------------------
  // fftw plans
  //--------------------
  bool plan_r2c_3d_created;
  T* plan_r2c_3d_in_pointer;
  fftw_plan plan_r2c_3d;

  bool plan_c2r_3d_created;
  T* plan_c2r_3d_out_pointer;
  fftw_plan plan_c2r_3d;

  // Grid size
  const int sizeX;
  const int sizeY;
  const int sizeZ;
  
  // Fourier transformed charge density
  CpuGrid<T> rhoK;

  // Green's function
  CpuGreensFunc<T>& greensFunc;
  
  void setupPlans(CpuGrid<T>& phi, const CpuGrid<T>& rho);

public:
  CpuFFTSolver(CpuGreensFunc<T>& greensFunc, const int sizeX, const int sizeY, const int sizeZ);
  ~CpuFFTSolver();
  void run(CpuGrid<T>& phi, const CpuGrid<T>& rho, const double sigma,
	   const double boxx, const double boxy, const double boxz);
};
#endif //CPUFFTSOLVER_H
