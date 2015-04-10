#include <iostream>
#include "TypesGSE.h"
#include "cuda/cuda_utils.h"
#include "cuda/XYZQ.h"
#include "cuda/Force.h"
#include "cuda/CudaPMERecip.h"

void testRandom_gpu(const int numCoord, const double L, const int ngrid, const double sigma, const int order,
		    xyzq_t<double>* xyzq, double& energy_SPME, double* forceXSPME, double* forceYSPME, double* forceZSPME) {

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
  
  double virial[9];
  energyVirial.copyToHost();
  cudaCheck(cudaDeviceSynchronize());
  energy_SPME = energyVirial.getEnergy("recip");
  //energy_self = energyVirial.getEnergy("self");
  energyVirial.getVirial(virial);

  std::cout << "energy_SPME = " << energy_SPME << std::endl;
  
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
