#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "CpuGSEr.h"
#ifdef USE_FFTW
#include "CpuGSEk.h"
#endif
#include "CpuGSEu.h"
#include "CpuEwaldRecip.h"
#include "CpuEwaldDirect.h"
#include "CpuGaussCharge.h"
#include "CpuGreensFuncG2.h"
#include "CpuLES.h"
#ifdef CUDA_ON
// CUDA
#include "cuda/cuda_utils.h"
#include "gpu_gse.h"
#endif

// Defines what we're computing:
#ifdef USE_FFTW
//#define USE_GSE_K
#endif

#ifdef CUDA_ON
#define USE_SPME
#endif

//#define USE_GSE_R
//#define USE_GSE_U
#define USE_LES

void test48K();
void testRandom(const int numCoord, const double L, const int ngrid, const double sigma, const int order, unsigned seed);
void testPair(const double r, const double L, const int ngrid, const double sigma);

int main(int argc, char *argv[]) {

  double r = 4.0;
  double L = 16.0;
  int M = -1;
  double sigma = 1.5;
  int order = 4;
  unsigned seed = 1234567890;
  
  bool arg_ok = true;
  int iarg = 1;
  while (iarg < argc) {
    if (strcmp(argv[iarg],"-r")==0) {
      iarg++;
      if (iarg == argc) {
	arg_ok = false;
	break;
      }
      sscanf(argv[iarg],"%lf",&r);
      iarg++;
    } else if (strcmp(argv[iarg],"-L")==0) {
      iarg++;
      if (iarg == argc) {
	arg_ok = false;
	break;
      }
      sscanf(argv[iarg],"%lf",&L);
      iarg++;
    } else if (strcmp(argv[iarg],"-M")==0) {
      iarg++;
      if (iarg == argc) {
	arg_ok = false;
	break;
      }
      sscanf(argv[iarg],"%d",&M);
      iarg++;
    } else if (strcmp(argv[iarg],"-sigma")==0) {
      iarg++;
      if (iarg == argc) {
	arg_ok = false;
	break;
      }
      sscanf(argv[iarg],"%lf",&sigma);
      iarg++;
    } else if (strcmp(argv[iarg],"-order")==0) {
      iarg++;
      if (iarg == argc) {
	arg_ok = false;
	break;
      }
      sscanf(argv[iarg],"%d",&order);
      iarg++;
    } else if (strcmp(argv[iarg],"-seed")==0) {
      iarg++;
      if (iarg == argc) {
	arg_ok = false;
	break;
      }
      sscanf(argv[iarg],"%u",&seed);
      iarg++;
    } else {
      std::cout << "Invalid input parameter " << argv[iarg] << std::endl;
      arg_ok = false;
      break;
    }
  }

  if (M < 0) M = (int)L;
  
  if (!arg_ok || r < 0.0 || L <= 0.0 || M <= 0 || (order != 4 && order != 6 && order != 8)) {
    std::cout << "Usage: cpu_gse -r r -L L -M M -sigma sigma -order order -seed seed"<< std::endl;
    return 1;
  }

#ifdef CUDA_ON
  std::vector<int> devices;
  start_gpu(1, 0, devices);
#endif
  
  //testRandom(100, L, M, sigma, order, seed);
  testPair(r, L, M, sigma);
  //test48K();

#ifdef CUDA_ON
  stop_gpu();
#endif
  
  return 0;
}

//
// Calculates root mean square error between two arrays
//
double Erms(const int n, const double* refx, const double* refy, const double* refz,
	    const double* x, const double* y, const double* z) {
  double rms = 0.0;
  double ref = 0.0;
  for (int i=0;i < n;i++) {
    double rx = refx[i];
    double ry = refy[i];
    double rz = refz[i];
    double dx = x[i] - rx;
    double dy = y[i] - ry;
    double dz = z[i] - rz;
    rms += dx*dx + dy*dy + dz*dz;
    ref += rx*rx + ry*ry + rz*rz;
  }
  return sqrt(rms/ref);
}

//
// Loads array data from file
//
template <typename T>
void loadArray(const int width, const int numLine, const char *filename, T *array) {
  std::ifstream file(filename);
  if (file.is_open()) {

    for (int i=0;i < numLine;i++) {
      for (int k=0;k < width;k++) {
	if (!(file >> array[i*width+k])) {
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

template<typename T>
void saveForce(const char *filename, const int numCoord, const T* forceX, const T* forceY, const T* forceZ) {
  std::ofstream file(filename);
  if (file.is_open()) {
    for (int i=0;i < numCoord;i++) {
      file << forceX[i] << " " << forceY[i] << " " << forceZ[i] << std::endl;
    }
  } else {
    std::cerr<<"Error opening file "<<filename<<std::endl;
    exit(1);
  }
}

template<typename T>
void readForce(const char *filename, const int numCoord, T* forceX, T* forceY, T* forceZ) {
  std::ifstream file(filename);
  if (file.is_open()) {
    for (int i=0;i < numCoord;i++) {
      file >> forceX[i] >> forceY[i] >> forceZ[i];
    }
  } else {
    std::cerr<<"Error opening file "<<filename<<std::endl;
    exit(1);
  }
}

template<typename T>
void saveMatrix(const char* filename,const T* L, const int sizeX, const int sizeY) {
  std::ofstream file(filename);
  if (file.is_open()) {
    for (int y=0;y < sizeY;y++) {
      for (int x=0;x < sizeX;x++) {
	file << L[x + y*sizeX] << " ";
      }
      file << std::endl;
    }
  } else {
    std::cerr << "saveMatrix<T>, Error opening file " << filename << std::endl;
    exit(1);
  }
}

//#######################################################################################
//#######################################################################################
//#######################################################################################
void test48K() {
  const int numCoord = 48191;
  const double sigma = 2.12132034355964;
  const double kappa = 0.4;
  const double lambdaSigma = 3.0;
  const double lambdaSigma1 = 3.0;
  const double boxx = 80.0;
  const double boxy = 80.0;
  const double boxz = 80.0;
  const int ngrid = 80;
  const int order = 4;
  
  xyzq_t<double>* xyzq = new xyzq_t<double>[numCoord];
  loadArray<double>(4, numCoord, "data/xyzq_48k.txt", (double *)xyzq);

#ifdef USE_GSE_K
  //--------------------------------------------------------------------------------------------------------
  // GSE-k
  //--------------------------------------------------------------------------------------------------------
  double* forceXGSEk = new double[numCoord];
  double* forceYGSEk = new double[numCoord];
  double* forceZGSEk = new double[numCoord];
  CpuGSEk<double,double> cpuGSEk(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEk.printInfo();
  cpuGSEk.spreadCharge1(numCoord, xyzq);
  cpuGSEk.solvePoisson();
  double energy_GSEk = cpuGSEk.calculateEnergy();
  std::cout << "energy_GSEk = " << energy_GSEk << std::endl;
  cpuGSEk.interpolateForce(numCoord, xyzq, forceXGSEk, forceYGSEk, forceZGSEk);
  saveForce<double>("forceGSEk.txt", numCoord, forceXGSEk, forceYGSEk, forceZGSEk);
#endif
  
#ifdef USE_GSE_U
  //--------------------------------------------------------------------------------------------------------
  // GSE-u
  //--------------------------------------------------------------------------------------------------------
  CpuGSEu<double,double> cpuGSEu(4, 1.1639, ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEu.spreadCharge1(numCoord, xyzq);
  cpuGSEu.solvePoisson();
  double energy_GSEu = cpuGSEu.calculateEnergy();
  std::cout << "energy_GSEu = " << energy_GSEu << std::endl;
#endif

#ifdef USE_LES
  //--------------------------------------------------------------------------------------------------------
  // LES
  //--------------------------------------------------------------------------------------------------------
  double* forceXLES = new double[numCoord];
  double* forceYLES = new double[numCoord];
  double* forceZLES = new double[numCoord];
  CpuLES<double,double> cpuLES(ngrid, ngrid, ngrid, sigma, boxx, boxy, boxz);
  //cpuLES.spreadCharge1(sigma/sqrt(2.0), 3.0, numCoord, xyzq, true);
  double gamma = 2.0;
  double sigmaA = sigma/sqrt(2.0)/sqrt(1+gamma*gamma);
  double sigmaB = gamma*sigmaA;
  printf("sigmaA = %lf sigmaB = %lf gamma = %lf\n",sigmaA,sigmaB,gamma);
  cpuLES.spreadCharge1(sigmaA, 3.0, numCoord, xyzq, true);
  cpuLES.spreadCharge2(sigmaB, 3.0);

  cpuLES.clearMagneticField();
  cpuLES.initElectricFieldJR();
  double err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  for (int i=0;i < 20000;i++) {
    cpuLES.integrate(2.0, 0.1, 1.0, 0.0);
  }
  double energy_LES = cpuLES.calcTotalEnergy();
  double energy_dip = cpuLES.calcDipoleEnergy();
  std::cout << "energy_LES = " << energy_LES-energy_dip  << " (" << energy_LES << " , " << energy_dip << ")" << std::endl;
  err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  printf("max(curlB) = %e\n",cpuLES.maxCurlB());
  printf("max(curlE) = %e\n",cpuLES.maxCurlE());

  cpuLES.chargeFluctuation(sigma/sqrt(2.0), 3.0, numCoord, xyzq, xyzq, true);
  cpuLES.interpolateForce(sigma/sqrt(2.0), 3.0, numCoord, xyzq, forceXLES, forceYLES, forceZLES);
  saveForce<double>("forceLES.txt", numCoord, forceXLES, forceYLES, forceZLES);
#endif

#ifdef USE_SPME
  //--------------------------------------------------------------------------------------------------------
  // SPME
  //--------------------------------------------------------------------------------------------------------
  double* forceXSPME = new double[numCoord];
  double* forceYSPME = new double[numCoord];
  double* forceZSPME = new double[numCoord];
  double* forceXSPMEdir = new double[numCoord];
  double* forceYSPMEdir = new double[numCoord];
  double* forceZSPMEdir = new double[numCoord];
  double energy_SPME;
  testRandom_gpu(numCoord, boxx, ngrid, sigma, order, xyzq, energy_SPME, forceXSPME, forceYSPME, forceZSPME,
		 NULL, NULL, NULL);
  saveForce<double>("forceSPME.txt", numCoord, forceXSPME, forceYSPME, forceZSPME);
  //saveForce<double>("forceSPMEdir.txt", numCoord, forceXSPMEdir, forceYSPMEdir, forceZSPMEdir);
#endif

  //--------------------------------------------------------------------------------------------------------  
  // Ewald
  //--------------------------------------------------------------------------------------------------------
  double* forceXew = new double[numCoord];
  double* forceYew = new double[numCoord];
  double* forceZew = new double[numCoord];
  double* forceXewdir = new double[numCoord];
  double* forceYewdir = new double[numCoord];
  double* forceZewdir = new double[numCoord];
  double energy_ew;
  //#define CALC_EWALD
#ifdef CALC_EWALD
  CpuEwaldRecip cpuEwald(sigma, 64, boxx, boxy, boxz);
  energy_ew = cpuEwald.calcForceEnergy(numCoord, xyzq, forceXew, forceYew, forceZew);
  saveForce<double>("forceEW64_48K.txt",numCoord, forceXew, forceYew, forceZew);
#else
  energy_ew = 1733.61;
  readForce<double>("forceEW64_48K.txt",numCoord, forceXew, forceYew, forceZew);
#endif
  std::cout << "energy_ew = " << energy_ew << std::endl;
  const int niblo14 = 48191;
  const int ninb14 = 61963;
  int *iblo14 = new int[niblo14];
  int *inb14 = new int[ninb14];
  loadArray<int>(1, niblo14, "data/iblo14_48k.txt", iblo14);
  loadArray<int>(1, ninb14, "data/inb14_48k.txt", inb14);
  CpuEwaldDirect cpuEWdirect(sigma, 9.0, boxx, boxy, boxz);
  cpuEWdirect.setupExclusions(niblo14, iblo14, inb14);
  double energy_ew_dir = cpuEWdirect.calcForceEnergy(numCoord, xyzq, forceXewdir, forceYewdir, forceZewdir);
  saveForce<double>("forceEWdir_48K.txt",numCoord, forceXewdir, forceYewdir, forceZewdir);
  delete [] iblo14;
  delete [] inb14;
  
#ifdef USE_GSE_K
  printf("rms(GSEk) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXGSEk, forceYGSEk, forceZGSEk),
	 fabs((energy_ew-energy_GSEk)/energy_ew));
#endif
#ifdef USE_LES
  printf("rms(LES) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXLES, forceYLES, forceZLES),
	 fabs((energy_ew-energy_LES)/energy_ew));
#endif
#ifdef USE_SPME
  printf("rms(SPME%d) = %e %e\n",order,Erms(numCoord, forceXew, forceYew, forceZew, forceXSPME, forceYSPME, forceZSPME),
	 fabs((energy_ew-energy_SPME)/energy_ew));
#endif

  // Add in direct contribution
  for (int i=0;i < numCoord;i++) {
#ifdef USE_LES
    forceXLES[i] += forceXewdir[i];
    forceYLES[i] += forceYewdir[i];
    forceZLES[i] += forceZewdir[i];
#endif
#ifdef USE_SPME
    forceXSPME[i] += forceXewdir[i];
    forceYSPME[i] += forceYewdir[i];
    forceZSPME[i] += forceZewdir[i];
#endif
    forceXew[i] += forceXewdir[i];
    forceYew[i] += forceYewdir[i];
    forceZew[i] += forceZewdir[i];
  }
#ifdef USE_LES
  energy_LES += energy_ew_dir;
#endif
#ifdef USE_SPME
  energy_SPME += energy_ew_dir;
#endif
  energy_ew += energy_ew_dir;
  
#ifdef USE_LES
  printf("rms(LES)  = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXLES, forceYLES, forceZLES),
	 fabs((energy_ew-(energy_LES-energy_dip))/energy_ew));
#endif
#ifdef USE_SPME
  printf("rms(SPME%d) = %e %e\n",order,Erms(numCoord, forceXew, forceYew, forceZew, forceXSPME, forceYSPME, forceZSPME),
	 fabs((energy_ew-energy_SPME)/energy_ew));
#endif

  delete [] xyzq;

#ifdef USE_GSE_K
  delete [] forceXGSEk;
  delete [] forceYGSEk;
  delete [] forceZGSEk;
#endif
  
#ifdef USE_LES
  delete [] forceXLES;
  delete [] forceYLES;
  delete [] forceZLES;
#endif

#ifdef USE_SPME
  delete [] forceXSPME;
  delete [] forceYSPME;
  delete [] forceZSPME;
  delete [] forceXSPMEdir;
  delete [] forceYSPMEdir;
  delete [] forceZSPMEdir;
#endif
  
  delete [] forceXew;
  delete [] forceYew;
  delete [] forceZew;

}

//#######################################################################################
//#######################################################################################
//#######################################################################################
void testRandom(const int numCoord, const double L, const int ngrid, const double sigma, const int order, unsigned seed) {
  const double kappa = 0.4;
  const double lambdaSigma = 3.0;
  const double lambdaSigma1 = 3.0;
  const double boxx = L;
  const double boxy = L;
  const double boxz = L;

  srand(seed);
  xyzq_t<double>* xyzq = new xyzq_t<double>[numCoord];
  for (int i=0;i < numCoord;i++) {
    xyzq[i].x = boxx*((double)rand()/((double)RAND_MAX + 1.0));
    xyzq[i].y = boxy*((double)rand()/((double)RAND_MAX + 1.0));
    xyzq[i].z = boxz*((double)rand()/((double)RAND_MAX + 1.0));
    xyzq[i].q = (double)((i%2)*2 - 1);
  }

#ifdef USE_GSE_K
  //--------------------------------------------------------------------------------------------------------
  // GSE-k
  //--------------------------------------------------------------------------------------------------------
  double* forceXGSEk = new double[numCoord];
  double* forceYGSEk = new double[numCoord];
  double* forceZGSEk = new double[numCoord];
  CpuGSEk<double,double> cpuGSEk(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEk.printInfo();
  cpuGSEk.spreadCharge1(numCoord, xyzq);
  cpuGSEk.solvePoisson();
  double energy_GSEk = cpuGSEk.calculateEnergy();
  std::cout << "energy_GSEk = " << energy_GSEk << std::endl;
  cpuGSEk.interpolateForce(numCoord, xyzq, forceXGSEk, forceYGSEk, forceZGSEk);
  saveForce<double>("forceGSEk.txt", numCoord, forceXGSEk, forceYGSEk, forceZGSEk);
#endif
  
#ifdef USE_GSE_R
  //--------------------------------------------------------------------------------------------------------
  // GSE-r
  //--------------------------------------------------------------------------------------------------------
  double* forceXGSEr = new double[numCoord];
  double* forceYGSEr = new double[numCoord];
  double* forceZGSEr = new double[numCoord];
  CpuGSEr<double,double> cpuGSEr(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEr.printInfo();
  cpuGSEr.spreadCharge1(numCoord, xyzq);
  cpuGSEr.spreadCharge2();
  cpuGSEr.solvePoisson();
  double energy_GSEr = cpuGSEr.calculateEnergy();
  std::cout << "energy_GSEr = " << energy_GSEr << std::endl;
  cpuGSEr.interpolateForce(numCoord, xyzq, forceXGSEr, forceYGSEr, forceZGSEr);
  saveForce<double>("forceGSEr.txt", numCoord, forceXGSEr, forceYGSEr, forceZGSEr);
#endif
  
#ifdef USE_GSE_U
  //--------------------------------------------------------------------------------------------------------
  // GSE-u
  //--------------------------------------------------------------------------------------------------------
  CpuGSEu<double,double> cpuGSEu(4, 2.0, ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEu.spreadCharge1(numCoord, xyzq);
  cpuGSEu.solvePoisson();
  double energy_GSEu = cpuGSEu.calculateEnergy();
  std::cout << "energy_GSEu = " << energy_GSEu << std::endl;
#endif
  
#ifdef USE_LES
  //--------------------------------------------------------------------------------------------------------
  // LES
  //--------------------------------------------------------------------------------------------------------
  double* forceXLES = new double[numCoord];
  double* forceYLES = new double[numCoord];
  double* forceZLES = new double[numCoord];
  CpuLES<double,double> cpuLES(ngrid, ngrid, ngrid, sigma, boxx, boxy, boxz);
    //#define USE_ONGRID_SPREAD
#ifdef USE_ONGRID_SPREAD
  double sigmaA = cpuGSEk.getSigma1();
  double sigmaB = cpuGSEk.getSigma2()/sqrt(2.0);
  cpuLES.spreadCharge1(sigmaA, 3.0, numCoord, xyzq);
#else
  cpuLES.spreadCharge1(sigma/sqrt(2.0), 3.0, numCoord, xyzq, true);
#endif
  //cpuLES.spreadCharge1(cpuGSEk.getSigma1(), cpuGSEk.getLambdaSigma1(), 2, xyzq);
  //cpuLES.spreadCharge2(cpuGSEk.getSigma2(), cpuGSEk.getLambdaSigma1());
  
  cpuLES.clearMagneticField();
  cpuLES.initElectricFieldJR();
  double err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  for (int i=0;i < 10000;i++) {
    cpuLES.integrate(2.0, 0.1, 1.0, 0.0);
  }
#ifdef USE_ONGRID_SPREAD
  double energy_LES = cpuLES.calcTotalEnergy(sigmaB, 3.0);
#else
  double energy_LES = cpuLES.calcTotalEnergy();
#endif
  double energy_dip = cpuLES.calcDipoleEnergy();
  std::cout << "energy_LES = " << energy_LES-energy_dip  << " (" << energy_LES << " , " << energy_dip << ")" << std::endl;
  err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  printf("max(curlB) = %e\n",cpuLES.maxCurlB());
  printf("max(curlE) = %e\n",cpuLES.maxCurlE());

#ifdef USE_ONGRID_SPREAD
  cpuLES.chargeFluctuation(sqrt(sigmaA*sigmaA + sigmaB*sigmaB), 3.0, numCoord, xyzq, xyzq, true);
  cpuLES.interpolateForce(sqrt(sigmaA*sigmaA + sigmaB*sigmaB), 3.0, numCoord, xyzq, forceXLES, forceYLES, forceZLES);
#else
  cpuLES.chargeFluctuation(sigma/sqrt(2.0), 3.0, numCoord, xyzq, xyzq, true);
  cpuLES.interpolateForce(sigma/sqrt(2.0), 3.0, numCoord, xyzq, forceXLES, forceYLES, forceZLES);
#endif
  
  saveForce<double>("forceLES.txt", numCoord, forceXLES, forceYLES, forceZLES);
#endif

#ifdef USE_SPME
  //--------------------------------------------------------------------------------------------------------
  // SPME
  //--------------------------------------------------------------------------------------------------------
  double* forceXSPME = new double[numCoord];
  double* forceYSPME = new double[numCoord];
  double* forceZSPME = new double[numCoord];
  double energy_SPME;
  testRandom_gpu(numCoord, L, ngrid, sigma, order, xyzq, energy_SPME, forceXSPME, forceYSPME, forceZSPME,
		 NULL, NULL, NULL);
  saveForce<double>("forceSPME.txt", numCoord, forceXSPME, forceYSPME, forceZSPME);
#endif
  
  //--------------------------------------------------------------------------------------------------------  
  // Ewald
  //--------------------------------------------------------------------------------------------------------
  double* forceXew = new double[numCoord];
  double* forceYew = new double[numCoord];
  double* forceZew = new double[numCoord];
  double* forceXewdir = new double[numCoord];
  double* forceYewdir = new double[numCoord];
  double* forceZewdir = new double[numCoord];
  CpuEwaldRecip cpuEWrecip(sigma, 12, boxx, boxy, boxz);
  double energy_ew = cpuEWrecip.calcForceEnergy(numCoord, xyzq, forceXew, forceYew, forceZew);
  CpuEwaldDirect cpuEWdirect(sigma, 9.0, boxx, boxy, boxz);
  double energy_ew_dir = cpuEWdirect.calcForceEnergy(numCoord, xyzq, forceXewdir, forceYewdir, forceZewdir);
  std::cout << "energy_ew = " << energy_ew << std::endl;
  saveForce<double>("forceEW.txt", numCoord, forceXew, forceYew, forceZew);
#ifdef USE_LES
  printf("rms(LES)  = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXLES, forceYLES, forceZLES),
	 fabs((energy_ew-(energy_LES-energy_dip))/energy_ew));
#endif
#ifdef USE_GSE_K
  printf("rms(GSEk) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXGSEk, forceYGSEk, forceZGSEk),
	 fabs((energy_ew-energy_GSEk)/energy_ew));
#endif
#ifdef USE_GSE_R
  printf("rms(GSEr) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXGSEr, forceYGSEr, forceZGSEr),
	 fabs((energy_ew-energy_GSEr)/energy_ew));
#endif
#ifdef USE_SPME
  printf("rms(SPME) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXSPME, forceYSPME, forceZSPME),
	 fabs((energy_ew-energy_SPME)/energy_ew));
#endif

  // Add in direct contribution
  for (int i=0;i < numCoord;i++) {
#ifdef USE_LES
    forceXLES[i] += forceXewdir[i];
    forceYLES[i] += forceYewdir[i];
    forceZLES[i] += forceZewdir[i];
#endif
    forceXew[i] += forceXewdir[i];
    forceYew[i] += forceYewdir[i];
    forceZew[i] += forceZewdir[i];
  }
#ifdef USE_LES
  energy_LES += energy_ew_dir;
#endif
  energy_ew += energy_ew_dir;
  
#ifdef USE_LES
  printf("rms(LES)  = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXLES, forceYLES, forceZLES),
	 fabs((energy_ew-(energy_LES-energy_dip))/energy_ew));
#endif

  delete [] forceXew;
  delete [] forceYew;
  delete [] forceZew;
  delete [] forceXewdir;
  delete [] forceYewdir;
  delete [] forceZewdir;
  
  delete [] xyzq;

#ifdef USE_SPME
  delete [] forceXSPME;
  delete [] forceYSPME;
  delete [] forceZSPME;
#endif
  
#ifdef USE_GSE_K
  delete [] forceXGSEk;
  delete [] forceYGSEk;
  delete [] forceZGSEk;
#endif

#ifdef USE_GSE_R
  delete [] forceXGSEr;
  delete [] forceYGSEr;
  delete [] forceZGSEr;
#endif

#ifdef USE_LES
  delete [] forceXLES;
  delete [] forceYLES;
  delete [] forceZLES;
#endif
}

//#######################################################################################
//#######################################################################################
//#######################################################################################
void testPair(const double r, const double L, const int ngrid, const double sigma) {
  const double kappa = 0.4;
  const double lambdaSigma = 3.0;
  const double lambdaSigma1 = 3.0;
  const double boxx = L;
  const double boxy = L;
  const double boxz = L;
  const int order = 4;
  
  xyzq_t<double> xyzq[2];
  double forceX[2];
  double forceY[2];
  double forceZ[2];

  /*
  xyzq[0].x = -r/2.0 + 0.5*boxx;
  xyzq[0].y = 0.5*boxy;
  xyzq[0].z = 0.5*boxz;
  xyzq[1].x = r/2.0 + 0.5*boxx;
  xyzq[1].y = 0.5*boxy;
  xyzq[1].z = 0.5*boxz;
  */
  
  /*
  xyzq[0].x = 0.5*boxx;
  xyzq[0].y = -r/2.0 + 0.5*boxy;
  xyzq[0].z = 0.5*boxz;
  xyzq[1].x = 0.5*boxx;
  xyzq[1].y = r/2.0 + 0.5*boxy;
  xyzq[1].z = 0.5*boxz;
  */

  /*
  xyzq[0].x = 0.5*boxx;
  xyzq[0].y = 0.5*boxy;
  xyzq[0].z = -r/2.0 + 0.5*boxz;
  xyzq[1].x = 0.5*boxx;
  xyzq[1].y = 0.5*boxy;
  xyzq[1].z = r/2.0 + 0.5*boxz;
  */

  /*
  double a = r/(2.0*sqrt(3.0));
  xyzq[0].x = -a + 0.5*L;
  xyzq[0].y = -a + 0.5*L;
  xyzq[0].z = -a + 0.5*L;
  xyzq[1].x = a + 0.5*L;
  xyzq[1].y = a + 0.5*L;
  xyzq[1].z = a + 0.5*L;
  */
  
  double a = r/sqrt(3.0);
  xyzq[0].x = 0.0;
  xyzq[0].y = 0.0;
  xyzq[0].z = 0.0;
  xyzq[1].x = a + 0.0;
  xyzq[1].y = a + 0.0;
  xyzq[1].z = a + 0.0;

  xyzq[0].q = -1.0;
  xyzq[1].q = 1.0;

#ifdef USE_GSE_R
  //--------------------------------------------------------------------------------------------------------
  // GSE-r
  //--------------------------------------------------------------------------------------------------------
  CpuGSEr<double,double> cpuGSEr(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEr.printInfo();
  cpuGSEr.spreadCharge1(2, xyzq);
  cpuGSEr.spreadCharge2();
  cpuGSEr.solvePoisson();
  double energy_GSEr = cpuGSEr.calculateEnergy();
  std::cout << "energy_GSEr " << energy_GSEr << std::endl;
  cpuGSEr.interpolateForce(2, xyzq, forceX, forceY, forceZ);
  printf("GSEr: %.12lf %.12lf %.12lf | %.12lf %.12lf %.12lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
#endif

#ifdef USE_GSE_K
  //--------------------------------------------------------------------------------------------------------
  // GSE-k
  //--------------------------------------------------------------------------------------------------------
  CpuGSEk<double,double> cpuGSEk(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEk.printInfo();
  cpuGSEk.spreadCharge1(2, xyzq);
  cpuGSEk.solvePoisson();
  double energy_GSEk = cpuGSEk.calculateEnergy();
  std::cout << "energy_GSEk " << energy_GSEk << std::endl;
  CpuGrid<double> Ex(ngrid, ngrid, ngrid);
  CpuGrid<double> Ey(ngrid, ngrid, ngrid);
  CpuGrid<double> Ez(ngrid, ngrid, ngrid);
  CpuGaussCharge<double,double> gaussCharge(ngrid, ngrid, ngrid);
  gaussCharge.calcElectricFieldOnGrid(boxx, boxy, boxz, cpuGSEk.getPhi(), Ex, Ey, Ez);
  Ex.save("ExGSEk.txt");
  Ey.save("EyGSEk.txt");
  Ez.save("EzGSEk.txt");
  cpuGSEk.getPhi().save("phi.txt");
  cpuGSEk.interpolateForce(2, xyzq, forceX, forceY, forceZ);
  printf("GSEk: %.12lf %.12lf %.12lf | %.12lf %.12lf %.12lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
#endif

#ifdef USE_GSE_U
  //--------------------------------------------------------------------------------------------------------
  // GSE-u
  //--------------------------------------------------------------------------------------------------------
  CpuGSEu<double,double> cpuGSEu(4, 1.168, ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEu.spreadCharge1(2, xyzq);
  cpuGSEu.solvePoisson();
  double energy_GSEu = cpuGSEu.calculateEnergy();
  std::cout << "energy_GSEu " << energy_GSEu << std::endl;
  cpuGSEu.interpolateForce(2, xyzq, forceX, forceY, forceZ);
  printf("GSEu: %.12lf %.12lf %.12lf | %.12lf %.12lf %.12lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
  cpuGSEu.getPhi().save("phiU.txt");
#endif
  
  //#define TEST_CHARGE_SPREAD
#ifdef TEST_CHARGE_SPREAD
  //--------------------------------------------------------------------------------------------------------
  // Test charge spreading
  // cpuGSEr.getRhoS()  = real-space final charge density
  // cpuGSEk.getRhoS()  = initial charge density
  // rhoS_k             = final charge density constructed using back-and-forth Fourier transform

  // Spread the charge directly using sigma = sqrt(sigma1^2 + sigma2^2) (= sigma_t)
  // This is used as the reference result
  CpuGrid<double> rhoS_ref(ngrid, ngrid, ngrid);
  CpuGaussCharge<double,double> gaussCharge(ngrid, ngrid, ngrid);
  double sigma_t = sqrt(cpuGSEk.getSigma1()*cpuGSEk.getSigma1() + cpuGSEk.getSigma2()*cpuGSEk.getSigma2());
  printf("sigma_t = %lf\n",sigma_t);
  gaussCharge.spreadChargeToGrid(sigma_t, 3.0*sqrt(2.0)*sigma_t, 2, xyzq,
				 boxx, boxy, boxz, rhoS_ref);
  double rhoS_ref_max = rhoS_ref.maxAbsValue();

  // Compare to GSE-r
  CpuGrid<double> res(ngrid, ngrid, ngrid);
#ifdef USE_GSE_R
  res.copy(cpuGSEr.getRhoS());
  res.sub(rhoS_ref);
  printf("RE|rhoS_r - rhoS_ref| = %e\n",res.maxAbsValue()/rhoS_ref_max);
#endif
  
  // Compare to GSE-k
  // NOTE: rhoS_k charge density is created by back-and-forth FFT of the correct Green's function
  CpuGrid<double> rhoS_k(ngrid, ngrid, ngrid);
  CpuGreensFunc<double>* greensFunc = new CpuGreensFuncG2<double>();
  CpuFFTSolver<double> fftSolver(*greensFunc, ngrid, ngrid, ngrid);
  fftSolver.run(rhoS_k, cpuGSEk.getRhoS(), cpuGSEk.getSigma2(), boxx, boxy, boxz);
  delete greensFunc;
  res.copy(rhoS_k);
  res.sub(rhoS_ref);
  printf("RE|rhoS_k - rhoS_ref| = %e\n",res.maxAbsValue()/rhoS_ref_max);

#ifdef USE_GSE_R
  cpuGSEr.getRhoS().save("rhoS_r.txt");
#endif
  rhoS_ref.save("rhoS_ref.txt");
  rhoS_k.save("rhoS_k.txt");
#endif

#ifdef USE_SPME
  //--------------------------------------------------------------------------------------------------------
  // SPME
  //--------------------------------------------------------------------------------------------------------
  double energy_SPME;
  testRandom_gpu(2, boxx, ngrid, sigma, order, xyzq, energy_SPME, forceX, forceY, forceZ, NULL, NULL, NULL);
  std::cout << "energy_SPME " << energy_SPME << std::endl;
  printf("SPME: %.12lf %.12lf %.12lf | %.12lf %.12lf %.12lf\n",
	 forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
#endif

#ifdef USE_LES
  //--------------------------------------------------------------------------------------------------------  
  // LES
  //--------------------------------------------------------------------------------------------------------  
  CpuLES<double,double> cpuLES(ngrid, ngrid, ngrid, sigma, boxx, boxy, boxz);
  //#define USE_ONGRID_SPREAD
#ifdef USE_ONGRID_SPREAD
  double sigmaA = cpuGSEk.getSigma1();
  double sigmaB = cpuGSEk.getSigma2()/sqrt(2.0);
  cpuLES.spreadCharge1(sigmaA, 3.0, 2, xyzq, true);
#else
  /*
  double gamma = 2.0;
  double sigmaA = sigma/sqrt(2.0)/sqrt(1+gamma*gamma);
  double sigmaB = gamma*sigmaA;
  printf("sigmaA = %lf sigmaB = %lf gamma = %lf\n",sigmaA,sigmaB,gamma);
  cpuLES.spreadCharge1(sigmaA, 3.0, 2, xyzq, true);
  cpuLES.spreadCharge2(sigmaB, 3.0);
  */
  cpuLES.spreadCharge1(sigma/sqrt(2.0), 2.5, 2, xyzq, true);
  cpuLES.getRho().save("rho.txt");
#endif
#ifdef TEST_CHARGE_SPREAD
  res.copy(cpuLES.getRho());
  res.sub(rhoS_ref);
  printf("RE|rhoS_LES - rhoS_ref| = %e\n",res.maxAbsValue()/rhoS_ref_max);
#endif
  cpuLES.clearMagneticField();
  
  //cpuLES.initElectricFieldJR();
  cpuLES.initElectricFieldLAPACK();
  
  cpuLES.getEx().save("ExLESinit.txt");
  cpuLES.getEy().save("EyLESinit.txt");
  cpuLES.getEz().save("EzLESinit.txt");
  double err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  for (int i=0;i < 10000;i++) {
      cpuLES.integrate(2.0, 0.1, 1.0, 0.0);
      //printf("Gauss, curlE, curlB = %e %e %e\n",cpuLES.checkGaussLaw(),cpuLES.maxCurlE(),cpuLES.maxCurlB());
  }
#ifdef USE_ONGRID_SPREAD
  double energy_LES = cpuLES.calcTotalEnergy(sigmaB, 3.0);
#else
  double energy_LES = cpuLES.calcTotalEnergy();
#endif
  double energy_dip = cpuLES.calcDipoleEnergy();
  std::cout << "energy_LES " << energy_LES-energy_dip  << " " << energy_LES << " " << energy_dip
	    << " " << cpuLES.calcDipoleEnergy(2, xyzq) << std::endl;
  err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  printf("max(curlB) = %e | energy_B = %e\n",cpuLES.maxCurlB(),cpuLES.calcTotalMagneticEnergy());
  printf("max(curlE) = %e\n",cpuLES.maxCurlE());

  cpuLES.getEx().save("ExLES.txt");
  cpuLES.getEy().save("EyLES.txt");
  cpuLES.getEz().save("EzLES.txt");

#ifdef USE_ONGRID_SPREAD
  cpuLES.chargeFluctuation(sqrt(sigmaA*sigmaA + sigmaB*sigmaB), 3.0, 2, xyzq, xyzq, true);
  cpuLES.interpolateForce(sqrt(sigmaA*sigmaA + sigmaB*sigmaB), 3.0, 2, xyzq, forceX, forceY, forceZ);
#else
  cpuLES.chargeFluctuation(sigma/sqrt(2.0), 2.5, 2, xyzq, xyzq, true);
  cpuLES.interpolateForce(sigma/sqrt(2.0), 2.5, 2, xyzq, forceX, forceY, forceZ);
#endif
  printf("LES(J): %.12lf %.12lf %.12lf | %.12lf %.12lf %.12lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);

#ifdef USE_GSE_R
  gaussCharge.calcElectricFieldOnGrid(boxx, boxy, boxz, cpuGSEr.getPhi(), Ex, Ey, Ez);
  Ex.save("ExGSEr.txt");
  Ey.save("EyGSEr.txt");
  Ez.save("EzGSEr.txt");
#endif

#endif
  //--------------------------------------------------------------------------------------------------------  
  // Ewald
  //--------------------------------------------------------------------------------------------------------
  double refX[2], refY[2], refZ[2];
  CpuEwaldRecip cpuEwald(sigma, 32*ngrid/boxx, boxx, boxy, boxz);
  double energy_ew = cpuEwald.calcForceEnergy(2, xyzq, refX, refY, refZ);
  std::cout << "energy_ew " << energy_ew << std::endl;
  printf("EWALD: %.12lf %.12lf %.12lf | %.12lf %.12lf %.12lf\n",refX[0],refY[0],refZ[0],refX[1],refY[1],refZ[1]);
#ifdef USE_LES
  printf("rms = %e\n",Erms(2, refX, refY, refZ, forceX, forceY, forceZ));
#endif
}

