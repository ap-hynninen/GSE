#include <cassert>
#include <cmath>
#include "CudaGaussCharge.h"
#include "cuda/gpu_utils.h"
//#include "cuda/cuda_utils.h"

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

//##########################################################################################
//##########################################################################################
//##########################################################################################

//
// Class creator
//
template <typename AT, typename CT>
CudaGaussCharge<AT, CT>::CudaGaussCharge(const int sizeX, const int sizeY, const int sizeZ) :
  sizeX(sizeX), sizeY(sizeY), sizeZ(sizeZ) {
}

//
// Class creator
//
template <typename AT, typename CT>
CudaGaussCharge<AT, CT>::~CudaGaussCharge() {
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

  int nthread = NTHREAD_PER_COORD*16;
  int nblock = (numCoord-1)/(nthread/NTHREAD_PER_COORD) + 1;
  int shmem_size = (nthread/NTHREAD_PER_COORD)*(2*nx+1 + 2*ny+1 + 2*nz+1)*sizeof(CT);
  printf("nthread=%d nblock=%d shmem_size=%d\n",nthread,nblock,shmem_size);
  spreadChargeAtomic<CT> <<< nblock, nthread, shmem_size >>>
    (numCoord, xyzq, sizeX, sizeY, sizeZ, hx, hy, hz,
     (CT)(1.0/hx), (CT)(1.0/hy), (CT)(1.0/hz), pref, one_sigma2sq,
     nx, ny, nz, rho.getDataPointer());

  cudaCheck(cudaGetLastError());  
}

//
// Explicit instances of CudaGaussCharge
//
template class CudaGaussCharge<float, float>;
template class CudaGaussCharge<double, float>;
template class CudaGaussCharge<double, double>;
