#include <iostream>
#include <vector>
#include "TypesGSE.h"
#include "cuda/cuda_utils.h"
#include "cuda/XYZQ.h"
#include "cuda/Force.h"
#include "cuda/CudaPMERecip.h"
#include "cuda/CudaNeighborList.h"
#include "cuda/CudaPMEDirectForce.h"

//
// Loads indices from file
//
template <typename T>
void load_ind(const int nind, const char *filename, const int n, T *ind) {
  std::ifstream file(filename);
  if (file.is_open()) {

    for (int i=0;i < n;i++) {
      for (int k=0;k < nind;k++) {
	if (!(file >> ind[i*nind+k])) {
	  std::cerr<<"Error reading file "<<filename<<std::endl;
	  exit(1);
	}
      }
    }

  } else {
    std::cerr<<"Error opening file "<<filename<<std::endl;
    exit(1);
  }

}

void testRandom_gpu(const int numCoord, const double L, const int ngrid, const double sigma, const int order,
		    xyzq_t<double>* xyzq, double& energy_SPME, double* forceXSPME, double* forceYSPME, double* forceZSPME,
		    double* forceXSPMEdir, double* forceYSPMEdir, double* forceZSPMEdir) {

  const FFTtype fft_type = BOX;
  float* forceXSPMEf = new float[numCoord];
  float* forceYSPMEf = new float[numCoord];
  float* forceZSPMEf = new float[numCoord];
  // Setup reciprocal vectors
  double recip[9];
  for (int i=0;i < 9;i++) recip[i] = 0;
  recip[0] = 1.0/L;
  recip[4] = 1.0/L;
  recip[8] = 1.0/L;

  float4* h_xyzq = new float4[numCoord];
  for (int i=0;i < numCoord;i++) {
    h_xyzq[i].x = xyzq[i].x;
    h_xyzq[i].y = xyzq[i].y;
    h_xyzq[i].z = xyzq[i].z;
    h_xyzq[i].w = xyzq[i].q;
  }
  
  CudaEnergyVirial energyVirial;
  
  XYZQ xyzq_gpu(numCoord);
  CudaPMERecip<int, float, float2> grid(ngrid, ngrid, ngrid, order, fft_type, 1, 0,
					energyVirial, "recip", "self");
  Force<float> force(numCoord);

  xyzq_gpu.set_xyzq(numCoord, h_xyzq);

  energyVirial.clear();
  
  grid.spread_charge(xyzq_gpu.xyzq, xyzq_gpu.ncoord, recip);
  grid.r2c_fft();
  grid.scalar_sum(recip, 1.0/(sigma*sqrt(2.0)), true, true);
  grid.c2r_fft();
  grid.gather_force(xyzq_gpu.xyzq, xyzq_gpu.ncoord, recip, force.stride(), force.xyz());
  force.getXYZ(forceXSPMEf, forceYSPMEf, forceZSPMEf);

    //---------------------------------------------------------------------------------------------
  // Direct
  //---------------------------------------------------------------------------------------------
  if (numCoord == 48191 && forceXSPMEdir != NULL && forceYSPMEdir != NULL && forceZSPMEdir != NULL) {
    float* forceXSPMEdirf = new float[numCoord];
    float* forceYSPMEdirf = new float[numCoord];
    float* forceZSPMEdirf = new float[numCoord];
  
    const int niblo14 = 48191;
    const int ninb14 = 61963;
    int *iblo14 = new int[niblo14];
    int *inb14 = new int[ninb14];
    load_ind<int>(1, "data/iblo14_48k.txt", niblo14, iblo14);
    load_ind<int>(1, "data/inb14_48k.txt", ninb14, inb14);
    CudaTopExcl topExcl(numCoord, iblo14, inb14);
    //CudaTopExcl topExcl(0, iblo14, inb14);
    
    // Create I vs. I interaction
    std::vector<int> numIntZone(8, 0);
    std::vector< std::vector<int> > intZones(8, std::vector<int>() );
    numIntZone.at(0) = 1;
    intZones.at(0).push_back(0);
    int zone_patom[9] = {0, 48191, 48191, 48191, 48191, 48191, 48191, 48191, 48191};
    
    XYZQ xyzq_sorted(numCoord, 32);
    int *loc2glo = NULL;
    allocate<int>(&loc2glo, numCoord);
    
    CudaNeighborList<32> nlist(topExcl, 1, 1, 1);
    nlist.registerList(numIntZone, intZones);
    //nlist.set_test(true);
    nlist.sort(0, zone_patom, xyzq_gpu.xyzq, xyzq_sorted.xyzq, loc2glo);
    nlist.build(0, zone_patom, L, L, L, 9.0, xyzq_sorted.xyzq, loc2glo);
    cudaCheck(cudaDeviceSynchronize());
    
    Force<long long int> force_fp(numCoord);
    force_fp.clear();
    
    CudaPMEDirectForce<long long int, float> dir(energyVirial, "vdw", "elec", "excl");
    dir.setup(L, L, L, 1.0/(sigma*sqrt(2.0)), 9.0, 8.0, 1.0, NONE, EWALD);
    dir.calc_force(xyzq_sorted.xyzq, nlist.getBuilder(0), false, false, force_fp.stride(), force_fp.xyz());
    force_fp.convert(force);
    cudaCheck(cudaDeviceSynchronize());
    
    delete [] iblo14;
    delete [] inb14;
    deallocate<int>(&loc2glo);

    force.getXYZ(forceXSPMEdirf, forceYSPMEdirf, forceZSPMEdirf);

    for (int i=0;i < numCoord;i++) {
      forceXSPMEdir[i] = forceXSPMEdirf[i];
      forceYSPMEdir[i] = forceYSPMEdirf[i];
      forceZSPMEdir[i] = forceZSPMEdirf[i];
    }

    delete [] forceXSPMEdirf;
    delete [] forceYSPMEdirf;
    delete [] forceZSPMEdirf;
  }
  //---------------------------------------------------------------------------------------------

  double virial[9];
  energyVirial.copyToHost();
  cudaCheck(cudaDeviceSynchronize());
  energy_SPME = energyVirial.getEnergy("recip");
  //energy_self = energyVirial.getEnergy("self");
  energyVirial.getVirial(virial);

  for (int i=0;i < numCoord;i++) {
    forceXSPMEf[i] = -forceXSPMEf[i];
    forceYSPMEf[i] = -forceYSPMEf[i];
    forceZSPMEf[i] = -forceZSPMEf[i];
  }

  delete [] h_xyzq;

  for (int i=0;i < numCoord;i++) {
    forceXSPME[i] = forceXSPMEf[i];
    forceYSPME[i] = forceYSPMEf[i];
    forceZSPME[i] = forceZSPMEf[i];
  }

  delete [] forceXSPMEf;
  delete [] forceYSPMEf;
  delete [] forceZSPMEf;

}
