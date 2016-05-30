#ifndef GPU_GSE_H
#define GPU_GSE_H

void testRandom_gpu(const int numCoord, const double L, const int ngrid, const double sigma, const int order,
		    xyzq_t<double>* xyzq, double& energy_SPME, double* forceXSPME, double* forceYSPME, double* forceZSPME,
		    double* forceXSPMEdir, double* forceYSPMEdir, double* forceZSPMEdir);

#endif //GPU_GSE_H
