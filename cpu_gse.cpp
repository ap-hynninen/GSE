#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "CpuGSEr.h"
#include "CpuGSEk.h"
#include "CpuGSEu.h"
#include "CpuEwaldRecip.h"
#include "CpuGaussCharge.h"
#include "CpuGreensFuncG2.h"
#include "CpuLES.h"

void test48K();
void testRandom(const int numCoord, const double L, const int ngrid, const double sigma, unsigned seed);
void testPair(const double r, const double L, const int ngrid, const double sigma);

int main(int argc, char *argv[]) {

  double r = 4.0;
  double L = 16.0;
  int M = -1;
  double sigma = 1.5;
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
  
  if (!arg_ok || r <= 0.0 || L <= 0.0 || M <= 0) {
    std::cout << "Usage: gpu_gse -r r -L L -M M -sigma sigma -seed seed"<< std::endl;
    return 1;
  }
  
  testRandom(100, L, M, sigma, seed);
  //testPair(r, L, M, sigma);
  //test48K();
  
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

void saveForce(const char *filename, const int numCoord, const double* forceX, const double* forceY, const double* forceZ) {
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

void test48K() {
  const int numCoord = 48191;
  const double sigma = 2.1213;
  const double kappa = 0.4;
  const double lambdaSigma = 3.0;
  const double lambdaSigma1 = 3.0;
  const double boxx = 80.0;
  const double boxy = 80.0;
  const double boxz = 80.0;
  const int ngrid = 80*2;
  
  xyzq_t<double>* xyzq = new xyzq_t<double>[numCoord];
  loadArray<double>(4, numCoord, "data/xyzq_48k.txt", (double *)xyzq);

  /*
  for (int i=0;i < numCoord;i++) {
    double x = xyzq[i].x;
    double y = xyzq[i].y;
    double z = xyzq[i].z;
    x += 0.5*boxx;
    y += 0.5*boxy;
    z += 0.5*boxz;
    if (x < 0.0) x += boxx;
    if (x >= boxx) x -= boxx;
    if (y < 0.0) y += boxy;
    if (y >= boxy) y -= boxy;
    if (z < 0.0) z += boxz;
    if (z >= boxz) z -= boxz;
    xyzq[i].x = x;
    xyzq[i].y = y;
    xyzq[i].z = z;
  }
  */
  
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
  saveForce("forceGSEk.txt", numCoord, forceXGSEk, forceYGSEk, forceZGSEk);

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

#define USE_LES
#ifdef USE_LES
  //--------------------------------------------------------------------------------------------------------
  // LES
  //--------------------------------------------------------------------------------------------------------
  double* forceXLES = new double[numCoord];
  double* forceYLES = new double[numCoord];
  double* forceZLES = new double[numCoord];
  CpuLES<double,double> cpuLES(ngrid, ngrid, ngrid, sigma, boxx, boxy, boxz);
  cpuLES.spreadCharge1(sigma/sqrt(2.0), 3.0, numCoord, xyzq);
  cpuLES.clearMagneticField();
  cpuLES.initElectricField();
  double err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  for (int i=0;i < 10000;i++) {
    if (i % 10 == 0) {
      printf("i=%d ",i);
      printf("max(curlB) = %e ",cpuLES.maxCurlB());
      printf("max(curlE) = %e | %e\n",cpuLES.maxCurlE(),cpuLES.calcTotalEnergy());
    }
    cpuLES.integrate(2.0, 0.1, 1.0);
  }
  double energy_LES = cpuLES.calcTotalEnergy();
  double energy_dip = cpuLES.calcDipoleEnergy(numCoord, xyzq);
  std::cout << "energy_LES = " << energy_LES-energy_dip  << " (" << energy_LES << " , " << energy_dip << ")" << std::endl;
  err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  printf("max(curlB) = %e\n",cpuLES.maxCurlB());
  printf("max(curlE) = %e\n",cpuLES.maxCurlE());

  cpuLES.chargeFluctuation(sigma/sqrt(2.0), 3.0, numCoord, xyzq, xyzq);
  cpuLES.interpolateForceJ(sigma/sqrt(2.0), 3.0, numCoord, xyzq, forceXLES, forceYLES, forceZLES);
  saveForce("forceLES.txt", numCoord, forceXLES, forceYLES, forceZLES);
#endif

  //--------------------------------------------------------------------------------------------------------  
  // Ewald
  //--------------------------------------------------------------------------------------------------------
  double* forceXew = new double[numCoord];
  double* forceYew = new double[numCoord];
  double* forceZew = new double[numCoord];
  CpuEwaldRecip cpuEwald(sigma, 24, boxx, boxy, boxz);
  double energy_ew = cpuEwald.calcForceEnergy(numCoord, xyzq, forceXew, forceYew, forceZew);
  std::cout << "energy_ew = " << energy_ew << std::endl;
  saveForce("forceEW.txt",numCoord, forceXew, forceYew, forceZew);

  printf("rms(GSEk) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXGSEk, forceYGSEk, forceZGSEk),
	 fabs((energy_ew-energy_GSEk)/energy_ew));
#ifdef USE_LES
  printf("rms(LES) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXLES, forceYLES, forceZLES),
	 fabs((energy_ew-energy_LES)/energy_ew));
#endif
  
  delete [] xyzq;

  delete [] forceXGSEk;
  delete [] forceYGSEk;
  delete [] forceZGSEk;

#ifdef USE_LES
  delete [] forceXLES;
  delete [] forceYLES;
  delete [] forceZLES;
#endif

  delete [] forceXew;
  delete [] forceYew;
  delete [] forceZew;

}

void testRandom(const int numCoord, const double L, const int ngrid, const double sigma, unsigned seed) {
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

#define USE_GSE_K
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
  saveForce("forceGSEk.txt", numCoord, forceXGSEk, forceYGSEk, forceZGSEk);
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
  
#define USE_LES
#ifdef USE_LES
  //--------------------------------------------------------------------------------------------------------
  // LES
  //--------------------------------------------------------------------------------------------------------
  double* forceXLES = new double[numCoord];
  double* forceYLES = new double[numCoord];
  double* forceZLES = new double[numCoord];
  CpuLES<double,double> cpuLES(ngrid, ngrid, ngrid, sigma, boxx, boxy, boxz);
  cpuLES.spreadCharge1(sigma/sqrt(2.0), lambdaSigma, numCoord, xyzq);
  cpuLES.clearMagneticField();
  cpuLES.initElectricFieldJR();
  double err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  for (int i=0;i < 10000;i++) {
    cpuLES.integrate(2.0, 0.1, 1.0);
  }
  double energy_LES = cpuLES.calcTotalEnergy();
  double energy_dip = cpuLES.calcDipoleEnergy(numCoord, xyzq);
  std::cout << "energy_LES = " << energy_LES-energy_dip  << " (" << energy_LES << " , " << energy_dip << ")" << std::endl;
  err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  printf("max(curlB) = %e\n",cpuLES.maxCurlB());
  printf("max(curlE) = %e\n",cpuLES.maxCurlE());

  cpuLES.chargeFluctuation(sigma/sqrt(2.0), lambdaSigma, numCoord, xyzq, xyzq);
  cpuLES.interpolateForceJ(sigma/sqrt(2.0), lambdaSigma, numCoord, xyzq, forceXLES, forceYLES, forceZLES);

  saveForce("forceLES.txt", numCoord, forceXLES, forceYLES, forceZLES);
#endif
  
  //--------------------------------------------------------------------------------------------------------  
  // Ewald
  //--------------------------------------------------------------------------------------------------------
  double* forceXew = new double[numCoord];
  double* forceYew = new double[numCoord];
  double* forceZew = new double[numCoord];
  CpuEwaldRecip cpuEwald(sigma, 12, boxx, boxy, boxz);
  double energy_ew = cpuEwald.calcForceEnergy(numCoord, xyzq, forceXew, forceYew, forceZew);
  std::cout << "energy_ew = " << energy_ew << std::endl;
  saveForce("forceEW.txt", numCoord, forceXew, forceYew, forceZew);
#ifdef USE_LES
  printf("rms(LES)  = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXLES, forceYLES, forceZLES),
	 fabs((energy_ew-energy_LES)/energy_ew));
#endif
#ifdef USE_GSE_K
  printf("rms(GSEk) = %e %e\n",Erms(numCoord, forceXew, forceYew, forceZew, forceXGSEk, forceYGSEk, forceZGSEk),
	 fabs((energy_ew-energy_GSEk)/energy_ew));
#endif

  delete [] forceXew;
  delete [] forceYew;
  delete [] forceZew;
  
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
}

void testPair(const double r, const double L, const int ngrid, const double sigma) {
  const double kappa = 0.4;
  const double lambdaSigma = 3.0;
  const double lambdaSigma1 = 3.0;
  const double boxx = L;
  const double boxy = L;
  const double boxz = L;
  
  xyzq_t<double> xyzq[2];
  double forceX[2];
  double forceY[2];
  double forceZ[2];

  xyzq[0].x = -r/2.0 + 0.5*boxx;
  xyzq[0].y = 0.5*boxy;
  xyzq[0].z = 0.5*boxz;
  xyzq[1].x = r/2.0 + 0.5*boxx;
  xyzq[1].y = 0.5*boxy;
  xyzq[1].z = 0.5*boxz;
  
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
  
  xyzq[0].q = -1.0;
  xyzq[1].q = 1.0;

#ifdef USE_GSE_R
  CpuGSEr<double,double> cpuGSEr(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEr.printInfo();
  cpuGSEr.spreadCharge1(2, xyzq);
  cpuGSEr.spreadCharge2();
  cpuGSEr.solvePoisson();
  double energy_GSEr = cpuGSEr.calculateEnergy();
  std::cout << "energy_GSEr = " << energy_GSEr << std::endl;
  cpuGSEr.interpolateForce(2, xyzq, forceX, forceY, forceZ);
  printf("GSEr: %lf %lf %lf | %lf %lf %lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
#endif
  
  CpuGSEk<double,double> cpuGSEk(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEk.printInfo();
  cpuGSEk.spreadCharge1(2, xyzq);
  cpuGSEk.solvePoisson();
  double energy_GSEk = cpuGSEk.calculateEnergy();
  std::cout << "energy_GSEk = " << energy_GSEk << std::endl;
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
  printf("GSEk: %lf %lf %lf | %lf %lf %lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);

  //#define USE_GSE_U
#ifdef USE_GSE_U
  //--------------------------------------------------------------------------------------------------------
  // GSE-u
  //--------------------------------------------------------------------------------------------------------
  CpuGSEu<double,double> cpuGSEu(4, 1.168, ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEu.spreadCharge1(2, xyzq);
  cpuGSEu.solvePoisson();
  double energy_GSEu = cpuGSEu.calculateEnergy();
  std::cout << "energy_GSEu = " << energy_GSEu << std::endl;
  cpuGSEu.interpolateForce(2, xyzq, forceX, forceY, forceZ);
  printf("GSEu: %lf %lf %lf | %lf %lf %lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
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

#define USE_LES
#ifdef USE_LES
  //--------------------------------------------------------------------------------------------------------  
  // LES
  //--------------------------------------------------------------------------------------------------------  
  CpuLES<double,double> cpuLES(ngrid, ngrid, ngrid, sigma, boxx, boxy, boxz);
  //cpuLES.spreadCharge1(cpuGSEk.getSigma1(), cpuGSEk.getLambdaSigma1(), 2, xyzq);
  //cpuLES.spreadCharge2(cpuGSEk.getSigma2(), cpuGSEk.getLambdaSigma1());
  cpuLES.spreadCharge1(sigma/sqrt(2.0), 3.0, 2, xyzq);
#ifdef TEST_CHARGE_SPREAD
  res.copy(cpuLES.getRho());
  res.sub(rhoS_ref);
  printf("RE|rhoS_LES - rhoS_ref| = %e\n",res.maxAbsValue()/rhoS_ref_max);
#endif
  cpuLES.clearMagneticField();
  cpuLES.initElectricField();
  double err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  for (int i=0;i < 1000;i++) {
    cpuLES.integrate(2.0, 0.1, 1.0);
  }
  double energy_LES = cpuLES.calcTotalEnergy();
  double energy_dip = cpuLES.calcDipoleEnergy(2, xyzq);
  std::cout << "energy_LES = " << energy_LES-energy_dip  << " (" << energy_LES << " , " << energy_dip << ")" << std::endl;
  err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  printf("max(curlB) = %e\n",cpuLES.maxCurlB());
  printf("max(curlE) = %e\n",cpuLES.maxCurlE());
 
  //cpuLES.interpolateForceVW(sigma/sqrt(2.0), 3.0, 2, xyzq, forceX, forceY, forceZ);
  //printf("LES(VW):  %lf %lf %lf | %lf %lf %lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
  //cpuLES.interpolateForceEF(sigma/sqrt(2.0), 3.0, 2, xyzq, forceX, forceY, forceZ);
  //cpuLES.interpolateForceEF(cpuGSEk.getSigma1(), cpuGSEk.getLambdaSigma1(), 2, xyzq, forceX, forceY, forceZ);
  printf("LES(EF): %lf %lf %lf | %lf %lf %lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
  cpuLES.getEx().save("ExLES.txt");
  cpuLES.getEy().save("EyLES.txt");
  cpuLES.getEz().save("EzLES.txt");

  cpuLES.chargeFluctuation(sigma/sqrt(2.0), 3.0, 2, xyzq, xyzq);
  cpuLES.interpolateForceJ(sigma/sqrt(2.0), 3.0, 2, xyzq, forceX, forceY, forceZ);
  printf("LES(J): %lf %lf %lf | %lf %lf %lf\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);

#ifdef USE_GSE_R
  CpuGrid<double> Ex(ngrid, ngrid, ngrid);
  CpuGrid<double> Ey(ngrid, ngrid, ngrid);
  CpuGrid<double> Ez(ngrid, ngrid, ngrid);
  gaussCharge.calcElectricFieldOnGrid(boxx, boxy, boxz, cpuGSEr.getPhi(), Ex, Ey, Ez);
  Ex.save("ExGSEr.txt");
  Ey.save("EyGSEr.txt");
  Ez.save("EzGSEr.txt");
#endif

  /*
  // Calculate electric field for GSEr from phi
  double ExPart[2];
  double EyPart[2];
  double EzPart[2];
#ifdef USE_GSE_R
  FILE* handle = fopen("E_GSEr.txt","wt");
  gaussCharge.interpolateElectricField(cpuGSEr.getSigma1(), cpuGSEr.getLambdaSigma1()*sqrt(2.0)*cpuGSEr.getSigma1(),
				       2, xyzq, boxx, boxy, boxz,
				       cpuGSEr.getPhi(), ExPart, EyPart, EzPart, handle);
  fclose(handle);
  printf("EPart(GSEr): %e %e %e | %e %e %e\n",ExPart[0],EyPart[0],EzPart[0],ExPart[1],EyPart[1],EzPart[1]);
#endif
  cpuLES.interpolateElectricField(cpuGSEk.getSigma1(), cpuGSEk.getLambdaSigma1(), 2, xyzq, ExPart, EyPart, EzPart);
  printf("EPart(LES):  %e %e %e | %e %e %e\n",ExPart[0],EyPart[0],EzPart[0],ExPart[1],EyPart[1],EzPart[1]);
  */
#endif
  //--------------------------------------------------------------------------------------------------------  
  // Ewald
  //--------------------------------------------------------------------------------------------------------
  double refX[2], refY[2], refZ[2];
  CpuEwaldRecip cpuEwald(sigma, 12*ngrid/boxx, boxx, boxy, boxz);
  double energy_ew = cpuEwald.calcForceEnergy(2, xyzq, refX, refY, refZ);
  std::cout << "energy_ew = " << energy_ew << std::endl;
  printf("EWALD: %lf %lf %lf | %lf %lf %lf\n",refX[0],refY[0],refZ[0],refX[1],refY[1],refZ[1]);
#ifdef USE_LES
  printf("rms = %e\n",Erms(2, refX, refY, refZ, forceX, forceY, forceZ));
#endif
}

/*
void testDHFR() {
  const int numCoord = 23558;
  const double sigma = 2.1213;
  const double kappa = 0.4;
  const double lambdaSigma = 3.0;
  const double lambdaSigma1 = 3.0;
  const double boxx = 62.23;
  const double boxy = 62.23;
  const double boxz = 62.23;
  const int ngrid = 64;
  
  xyzq_t<double>* xyzq = new xyzq_t<double>[numCoord];
  loadArray<double>(4, numCoord, "data/xyzq.txt", (double *)xyzq);
  double* forceX = new double[numCoord];
  double* forceY = new double[numCoord];
  double* forceZ = new double[numCoord];
  
  // Fix up the system
  double sum_q = 0.0;
  for (int i=0;i < numCoord;i++) sum_q += xyzq[i].q;
  sum_q /= (double)numCoord;
  printf("sum_q = %e\n",sum_q);
  for (int i=0;i < numCoord;i++) xyzq[i].q -= sum_q;  
  sum_q = 0.0;
  for (int i=0;i < numCoord;i++) sum_q += xyzq[i].q;
  printf("sum_q = %e\n",sum_q);
  
  CpuGSEk<double,double> cpuGSEk(ngrid, ngrid, ngrid, sigma, kappa, lambdaSigma, lambdaSigma1, boxx, boxy, boxz);
  cpuGSEk.printInfo();
  cpuGSEk.spreadCharge1(numCoord, xyzq);
  cpuGSEk.solvePoisson();
  double energy_GSEk = cpuGSEk.calculateEnergy();
  std::cout << "energy_GSEk = " << energy_GSEk << std::endl;

  CpuLES<double,double> cpuLES(ngrid, ngrid, ngrid, sigma, boxx, boxy, boxz);
  cpuLES.spreadCharge1(sigma/sqrt(2.0), 3.0, numCoord, xyzq);
  cpuLES.clearMagneticField();
  cpuLES.initElectricField();
  double err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  for (int i=0;i < 10000;i++) {
    cpuLES.integrate(2.0, 0.1, 1.0);
  }
  double energy_LES = cpuLES.calcTotalEnergy();
  double energy_dip = cpuLES.calcDipoleEnergy(numCoord, xyzq);
  std::cout << "energy_LES = " << energy_LES-energy_dip  << " (" << energy_LES << " , " << energy_dip << ")" << std::endl;
  err = cpuLES.checkGaussLaw();
  printf("Error in Gauss Law = %e\n",err);
  printf("max(curlB) = %e\n",cpuLES.maxCurlB());
  printf("max(curlE) = %e\n",cpuLES.maxCurlE());

  cpuGSEk.interpolateForce(numCoord, xyzq, forceX, forceY, forceZ);
  printf("GSEk: %e %e %e | %e %e %e\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
  
  cpuLES.interpolateForceVW(sigma/sqrt(2.0), 3.0, numCoord, xyzq, forceX, forceY, forceZ);
  printf("LES:  %e %e %e | %e %e %e\n",forceX[0],forceY[0],forceZ[0],forceX[1],forceY[1],forceZ[1]);
  
  delete [] xyzq;
  delete [] forceX;
  delete [] forceY;
  delete [] forceZ;
}
*/
