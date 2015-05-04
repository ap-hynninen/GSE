#include <cassert>
#include <cmath>
#include "CudaGaussCharge.h"
#include "cuda/gpu_utils.h"
#include "cuda/cuda_utils.h"

#define USE_ATOMIC

//
// Simple atomic version
//
#define NTHREAD_PER_COORD 32
template <typename CT>
__global__ void spreadChargeAtomic(const int numCoord,
				   const xyzq_t<CT> *xyzq,
				   const int sizeX, const int sizeY, const int sizeZ,
				   const CT hx, const CT hy, const CT hz,
				   const CT inv_hx, const CT inv_hy, const CT inv_hz,
				   const CT pref, const CT one_sigma2sq,
				   const int nx, const int ny, const int nz,
				   CT* rho) {
  // Shared memory requirement: blockDim.x/NTHREAD_PER_COORD*(2*nx+1 + 2*ny+1 + 2*nz+1)*sizeof(CT)
  extern __shared__ char sh_buf[];

  // Setup shared memory pointers
  int sh_ind = threadIdx.x/NTHREAD_PER_COORD*(2*nx+1 + 2*ny+1 + 2*nz+1)*sizeof(CT);
  CT* sh_wx = (CT *)&sh_buf[sh_ind];
  sh_ind += (2*nx+1)*sizeof(CT);
  CT* sh_wy = (CT *)&sh_buf[sh_ind];
  sh_ind += (2*ny+1)*sizeof(CT);
  CT* sh_wz = (CT *)&sh_buf[sh_ind];

  const int icoord = (threadIdx.x + blockIdx.x*blockDim.x)/NTHREAD_PER_COORD;
  if (icoord < numCoord) {
    xyzq_t<CT> xyzqi = xyzq[icoord];
    int ixc = rintf(xyzqi.x*inv_hx) - nx;
    int iyc = rintf(xyzqi.y*inv_hy) - ny;
    int izc = rintf(xyzqi.z*inv_hz) - nz;
    const CT qi = xyzqi.q*pref;
    const CT dx = ixc*hx - xyzqi.x;
    const CT dy = iyc*hy - xyzqi.y;
    const CT dz = izc*hz - xyzqi.z;
    // Compute 1D weights
    for (int i=threadIdx.x % NTHREAD_PER_COORD;i <= 2*nx;i+=NTHREAD_PER_COORD) {
      CT x = dx + i*hx;
      sh_wx[i] = expf(-x*x*one_sigma2sq);
    }
    for (int i=threadIdx.x % NTHREAD_PER_COORD;i <= 2*ny;i+=NTHREAD_PER_COORD) {
      CT y = dy + i*hy;
      sh_wy[i] = expf(-y*y*one_sigma2sq);
    }
    for (int i=threadIdx.x % NTHREAD_PER_COORD;i <= 2*nz;i+=NTHREAD_PER_COORD) {
      CT z = dz + i*hz;
      sh_wz[i] = qi*expf(-z*z*one_sigma2sq);
    }
    // --- Warp level synchronization here ---
    const int nxw = (2*nx+1);
    const int nxyw = (2*nx+1)*(2*ny+1);
    const int nxyzw = (2*nx+1)*(2*ny+1)*(2*nz+1);
    // Make sure (ixc, iyc, izc) are non-negative
    ixc = (ixc + sizeX) % sizeX;
    iyc = (iyc + sizeY) % sizeY;
    izc = (izc + sizeZ) % sizeZ;
    for (int id=threadIdx.x % NTHREAD_PER_COORD;id < nxyzw;id+=NTHREAD_PER_COORD) {
      // Compute index to the weight table (i,j,k):
      // i=0...2*nx, j=0...2*ny, k=0...2*nz
      int t = id;
      int k = t/nxyw;
      t -= k*nxyw;
      int j = t/nxw;
      int i = t - j*nxw;

      // Compute weighted charge
      CT qw = sh_wx[i]*sh_wy[j]*sh_wz[k];

      // Compute index to the rho table
      int ix = (ixc + i) % sizeX;
      int iy = (iyc + j) % sizeY;
      int iz = (izc + k) % sizeZ;
      const int p = ix + (iy + iz*sizeY)*sizeX;
      atomicAdd(&rho[p], qw);
    }
  }
  
}

//
// Map coordinates to a 3d grid
//
template <typename CT>
__global__ void mapCoordToGrid(const int numCoord, const xyzq_t<CT> *xyzq,
			       const int sizeX, const int sizeY, const int sizeZ,
			       const CT hx, const CT hy, const CT hz,
			       const CT inv_hx, const CT inv_hy, const CT inv_hz,
			       xyzq_t<CT> *xyzqGrid, xyzq_t<CT> *xyzqOverflow,
			       int* pOverflow, int* coordP) {
  const int icoord = threadIdx.x + blockIdx.x*blockDim.x;
  if (icoord < numCoord) {
    xyzq_t<CT> xyzqi = xyzq[icoord];
    int ixc = rintf(xyzqi.x*inv_hx);// - nx;
    int iyc = rintf(xyzqi.y*inv_hy);// - ny;
    int izc = rintf(xyzqi.z*inv_hz);// - nz;    
    //xyzqi.x = ixc*hx - xyzqi.x;
    //xyzqi.y = iyc*hy - xyzqi.y;
    //xyzqi.z = izc*hz - xyzqi.z;
    // Make sure (ixc, iyc, izc) are non-negative
    ixc = (ixc + 8*sizeX) % sizeX;
    iyc = (iyc + 8*sizeY) % sizeY;
    izc = (izc + 8*sizeZ) % sizeZ;
    const int p = ixc + (iyc + izc*sizeY)*sizeX;
    coordP[icoord] = p;
    if (atomicCAS(&xyzqGrid[p], 0, __float_as_int(xyzqi.q)) != 0) {
      xyzqGrid[p] = xyzqi;
    } else {
      xyzqOverflow[atomicAdd(pOverflow, 1)] = xyzqi;
    }
  }
}

//
// Accumulate charges
//
template <typename CT>
__global__ void accumCharge(const xyzq_t<CT>* xyzqGrid, const int* coordP,
			    const int sizeX, const int sizeY, const int sizeZ,
			    CT* rho) {
  // Shared memory requirement: sizeof(CT)
  extern __shared__ char sh_buf[];
  CT* sh_w = &sh_buf[0];

  sh_w[threadIdx.x] = (CT)0;
  
  const int nyw = 2*ny+1;
  const int nyzw = (2*ny+1)*(2*nz+1);
  for (int ix=threadIdx.x;ix < sizeX;ix+=blockDim.x) {
    for (int id=0;id < nxyw;id++) {
      int iz = id/nyw;
      int iy = id-iz*nyw;
      iz = blockIdx.z - nz;
      iy = blockIdx.y - ny;
      int p = ix + (iy + iz*sizeY)*sizeX;
      xyzq_t<CT> xyzqi = xyzqGrid[p];
      if (xyzqi.q != (CT)0) {
	int ixc = rintf(xyzqi.x*inv_hx) - nx;
	int iyc = rintf(xyzqi.y*inv_hy) - ny;
	int izc = rintf(xyzqi.z*inv_hz) - nz;    
	CT dx = ixc*hx - xyzqi.x;
	CT dy = iyc*hy - xyzqi.y;
	CT dz = izc*hz - xyzqi.z;
	CT r2 = dx*dx + dy*dy + dz*dz;
	CT q_pref = xyzqi.q*pref;
	// Make sure (ixc, iyc, izc) are non-negative
	ixc = (ixc + 8*sizeX) % sizeX;
	iyc = (iyc + 8*sizeY) % sizeY;
	izc = (izc + 8*sizeZ) % sizeZ;
	sh_w[threadIdx.x] += q_pref*expf(-r2*one_sigma2sq);
      }
    }
  }
  
}

//
// Add in overflow charges
//
template <typename CT>
__global__ void overflowCharge(const xyzq_t<CT>* xyzqOverflow, int* pOverflow) {
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  if (i < *pOverflow) {
    xyzq_t<CT> xyzqi = xyzqOverflow[i];
    int ixc = rintf(xyzqi.x*inv_hx) - nx;
    int iyc = rintf(xyzqi.y*inv_hy) - ny;
    int izc = rintf(xyzqi.z*inv_hz) - nz;
    const CT qi = xyzqi.q*pref;
    const CT dx = ixc*hx - xyzqi.x;
    const CT dy = iyc*hy - xyzqi.y;
    const CT dz = izc*hz - xyzqi.z;
    // Compute 1D weights
    for (int i=threadIdx.x % NTHREAD_PER_COORD;i <= 2*nx;i+=NTHREAD_PER_COORD) {
      CT x = dx + i*hx;
      sh_wx[i] = expf(-x*x*one_sigma2sq);
    }
    for (int i=threadIdx.x % NTHREAD_PER_COORD;i <= 2*ny;i+=NTHREAD_PER_COORD) {
      CT y = dy + i*hy;
      sh_wy[i] = expf(-y*y*one_sigma2sq);
    }
    for (int i=threadIdx.x % NTHREAD_PER_COORD;i <= 2*nz;i+=NTHREAD_PER_COORD) {
      CT z = dz + i*hz;
      sh_wz[i] = qi*expf(-z*z*one_sigma2sq);
    }
    // --- Warp level synchronization here ---
    const int nxw = (2*nx+1);
    const int nxyw = (2*nx+1)*(2*ny+1);
    const int nxyzw = (2*nx+1)*(2*ny+1)*(2*nz+1);
    // Make sure (ixc, iyc, izc) are non-negative
    ixc = (ixc + sizeX) % sizeX;
    iyc = (iyc + sizeY) % sizeY;
    izc = (izc + sizeZ) % sizeZ;
    for (int id=threadIdx.x % NTHREAD_PER_COORD;id < nxyzw;id+=NTHREAD_PER_COORD) {
      // Compute index to the weight table (i,j,k):
      // i=0...2*nx, j=0...2*ny, k=0...2*nz
      int t = id;
      int k = t/nxyw;
      t -= k*nxyw;
      int j = t/nxw;
      int i = t - j*nxw;

      // Compute weighted charge
      CT qw = sh_wx[i]*sh_wy[j]*sh_wz[k];

      // Compute index to the rho table
      int ix = (ixc + i) % sizeX;
      int iy = (iyc + j) % sizeY;
      int iz = (izc + k) % sizeZ;
      const int p = ix + (iy + iz*sizeY)*sizeX;
      atomicAdd(&rho[p], qw);
    }
  }
}

//##########################################################################################
//##########################################################################################
//##########################################################################################

//
// Class creator
//
template <typename AT, typename CT>
CudaGaussCharge<AT, CT>::CudaGaussCharge(const int sizeX, const int sizeY, const int sizeZ) :
  sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ) {
#ifndef USE_ATOMIC
  allocate< xyzq_t<CT> >(&xyzqGrid, sizeX*sizeY*sizeZ);
  xyzqOverflowLen = pow(sizeX*sizeY*sizeZ,1.0/3.0) + 1000;
  allocate< xyzq_t<CT> >(&xyzqOverflow, xyzqOverflowLen);
  allocate<int>(&pOverflow, 1);
#endif
}

//
// Class creator
//
template <typename AT, typename CT>
CudaGaussCharge<AT, CT>::~CudaGaussCharge() {
#ifndef USE_ATOMIC
  deallocate< xyzq_t<CT> >(&xyzqGrid);
  deallocate< xyzq_t<CT> >(&xyzqOverflow);
  deallocate< int >(&pOverflow);
#endif
}

//
// Particle to grid charge spreading
//
template <typename AT, typename CT>
void CudaGaussCharge<AT, CT>::spreadChargeToGrid(const double sigma, const double rcut,
						 const int numCoord, const xyzq_t<CT> *xyzq,
						 const double boxx, const double boxy, const double boxz,
						 CudaGrid<CT>& rho) {

    // Sanity checks
  assert(sizeX == rho.getSizeX());
  assert(sizeY == rho.getSizeY());
  assert(sizeZ == rho.getSizeZ());
  
  const CT hx = (CT)(boxx/(double)sizeX);
  const CT hy = (CT)(boxy/(double)sizeY);
  const CT hz = (CT)(boxz/(double)sizeZ);
  const int nx = (int)ceil(rcut/hx);
  const int ny = (int)ceil(rcut/hy);
  const int nz = (int)ceil(rcut/hz);

  printf("CudaGaussCharge::spreadChargeToGrid, rcut=%lf nx,ny,nz=%d %d %d\n",rcut,nx,ny,nz);

  // Clear grid
  rho.clear();
  
  CT one_sigma2sq = (CT)(1.0/(2.0*sigma*sigma));
  CT pref = (CT)pow(2.0*pi_dbl*sigma*sigma,-3.0/2.0);

#ifdef USE_ATOMIC
  int nthread = NTHREAD_PER_COORD*16;
  int nblock = (numCoord-1)/(nthread/NTHREAD_PER_COORD) + 1;
  int shmem_size = (nthread/NTHREAD_PER_COORD)*(2*nx+1 + 2*ny+1 + 2*nz+1)*sizeof(CT);
  printf("nthread=%d nblock=%d shmem_size=%d\n",nthread,nblock,shmem_size);
  spreadChargeAtomic<CT> <<< nblock, nthread, shmem_size >>>
    (numCoord, xyzq, sizeX, sizeY, sizeZ, hx, hy, hz,
     (CT)(1.0/hx), (CT)(1.0/hy), (CT)(1.0/hz), pref, one_sigma2sq,
     nx, ny, nz, rho.getDataPointer());
  cudaCheck(cudaGetLastError());  
#else
  clear_gpu_array< xyzq_t<CT> >(xyzqGrid, sizeX*sizeY*sizeZ);
  clear_gpu_array<int>(pOverflow, 1);
  int nthread = 512;
  int nblock = (numCoord-1)/nthread + 1;
  mapCoordToGrid<<< nblock, nthread >>>
    (numCoord, xyzq, sizeX, sizeY, sizeZ, hx, hy, hz,
     (CT)(1.0/hx), (CT)(1.0/hy), (CT)(1.0/hz), xyzqGrid, xyzqOverflow, pOverflow);
  cudaCheck(cudaGetLastError());  
#endif

}

//
// Explicit instances of CudaGaussCharge
//
template class CudaGaussCharge<float, float>;
template class CudaGaussCharge<double, float>;
template class CudaGaussCharge<double, double>;
