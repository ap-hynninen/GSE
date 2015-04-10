#ifndef CPUMULTIGRIDSOLVER_H
#define CPUMULTIGRIDSOLVER_H
//
// Multigrid solver
// (c) Antti-Pekka Hynninen, Feb 2015
// aphynninen@hotmail.com
//

#include <vector>
#include "CpuGrid.h"

template <typename T>
class CpuMultiGridSolver {
private:

  // Laplace representations
  enum {SEVEN_POINT, HERMITIAN};
  
  // Number of levels
  int numLevel;
  
  // Number of pre and post relaxation sweeps
  int numPreRelax;
  int numPostRelax;

  //-------------------------------------
  // Grids, finest first
  //-------------------------------------
  // rho
  std::vector< CpuGrid<T>* > rhoVec;
  // phi
  std::vector< CpuGrid<T>* > phiVec;
  // residual
  std::vector< CpuGrid<T>* > resVec;

  void relaxRBGS(CpuGrid<T>& phi, const CpuGrid<T>& rho,
		 const double boxx, const double boxy, const double boxz,
		 const int laplaceRep);
  void relaxGS(CpuGrid<T>& phi, const CpuGrid<T>& rho,
	       const double boxx, const double boxy, const double boxz,
	       const int laplaceRep);
  void solveSmall(CpuGrid<T>& phi, const CpuGrid<T>& rho);
  void fineToCoarse(const CpuGrid<T>& rhoFine, CpuGrid<T>& rhoCoarse);
  void coarseToFine(const CpuGrid<T>& rhoCoarse, CpuGrid<T>& rhoFine);
  
public:
  CpuMultiGridSolver(const int sizeFinestX, const int sizeFinestY, const int sizeFinestZ);
  ~CpuMultiGridSolver();
  
  void run(CpuGrid<T>& phi, const CpuGrid<T>& rho,
	   const double boxx, const double boxy, const double boxz);

  void calculateResidual(const CpuGrid<T>& phi, const CpuGrid<T>& rho, CpuGrid<T>& res,
			 const double boxx, const double boxy, const double boxz,
			 const int laplaceRep);

  void buildLaplace(const CpuGrid<T>& phi, T* L, const double boxx, const double boxy, const double boxz);
  
};

#endif // CPUMULTIGRIDSOLVER_H
