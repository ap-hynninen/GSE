
# Detect OS
ifeq ($(shell uname -a|grep Linux|wc -l), 1)
OS = linux
endif

ifeq ($(shell uname -a|grep titan|wc -l), 1)
OS = titan
endif

ifeq ($(shell uname -a|grep Darwin|wc -l), 1)
OS = osx
endif

YES := $(shell which make | wc -l 2> /dev/null)
NO := $(shell which pikaboo | wc -l 2> /dev/null)

# Set optimization level
#OPTLEV = -g
OPTLEV = -O3

# Detect CUDA, Intel compiler, and MPI
ifeq ($(OS),titan)
MPI_FOUND := $(YES)
else
MPI_FOUND := $(shell which mpicc | wc -l 2> /dev/null)
endif

CUDA_COMPILER := $(shell which nvcc | wc -l 2> /dev/null)
INTEL_COMPILER := $(shell which icc | wc -l 2> /dev/null)

# Detect FFTW
FFTW_FOUND := $(shell which fftw-wisdom | wc -l 2> /dev/null)
ifeq ($(FFTW_FOUND), $(YES))
FFTWROOT = $(subst /bin/,,$(dir $(shell which fftw-wisdom)))
endif

#---------------------------------------------------------------
ifeq ($(MPI_FOUND), $(YES))

DEFS = -D USE_MPI

ifeq ($(OS),titan)
CC = CC
CL = CC
else  # ifeq ($(OS),titan)

ifeq ($(INTEL_COMPILER), $(YES))
CC = mpicc
CL = mpicc
DEFS += -D MPICH_IGNORE_CXX_SEEK
else
CC = mpic++
CL = mpic++
endif

endif  # ifeq ($(OS),titan)

else   # MPI_FOUND

DEFS = -D DONT_USE_MPI

ifeq ($(INTEL_COMPILER), $(YES))
CC = icc
CL = icc
DEFS += -D MPICH_IGNORE_CXX_SEEK
else
CC = g++
CL = g++
endif

endif  # MPI_FOUND

# NOTE: CUDA texture objects require Kepler GPU + CUDA 5.0
DEFS += -D USE_TEXTURE_OBJECTS

OBJS_GSE = cpu_gse.o CpuGSE.o CpuGSEr.o CpuGSEk.o CpuMultiGridSolver.o \
           CpuFFTSolver.o CpuEwaldRecip.o CpuGaussCharge.o CpuGreensFuncGSE.o \
	   CpuGreensFuncG2.o CpuLES.o CpuGSEu.o

OBJS_GPU_CONV = gpu_conv.o cuda_utils.o CudaConvolution.o CpuConvolution.o

ifeq ($(MPI_FOUND), $(YES))
OBJS = $(OBJS_GSE)
endif
ifeq ($(CUDA_COMPILER), $(YES))
OBJS += $(OBJS_GPU_CONV)
endif

ifeq ($(CUDA_COMPILER), $(YES))
CUDAROOT = $(subst /bin/,,$(dir $(shell which nvcc)))
endif

ifeq ($(MPI_FOUND), $(YES))
ifeq ($(OS),titan)
# NOTE: Assumes we're using Intel compiler
#MPIROOT = $(subst /bin/,,$(dir $(shell which icc)))
MPIROOT = /opt/cray/mpt/6.3.0/gni/mpich2-intel/130
else
MPIROOT = $(subst /bin/,,$(dir $(shell which mpicc)))
endif
endif

ifeq ($(INTEL_COMPILER), $(YES))
OPENMP_OPT = -openmp
else
OPENMP_OPT = -fopenmp
endif

ifeq ($(CUDA_COMPILER), $(YES))
GENCODE_SM20  := -gencode arch=compute_20,code=sm_20
GENCODE_SM30  := -gencode arch=compute_30,code=sm_30
GENCODE_SM35  := -gencode arch=compute_35,code=sm_35
GENCODE_SM50  := -gencode arch=compute_50,code=sm_50
GENCODE_FLAGS := $(GENCODE_SM20) $(GENCODE_SM30) $(GENCODE_SM35)
# See if CUDA compiler supports compute 5.0
ifneq ($(shell nvcc --help|grep compute_50|wc -l), 0)
GENCODE_FLAGS += $(GENCODE_SM50)
endif
endif

# CUDA_CFLAGS = flags for compiling CUDA API calls using c compiler
# NVCC_CFLAGS = flags for nvcc compiler
# CUDA_LFLAGS = flags for linking with CUDA

CUDA_CFLAGS = -I${CUDAROOT}/include $(OPTLEV) $(OPENMP_OPT) -std=c++0x
NVCC_CFLAGS = $(OPTLEV) -lineinfo -fmad=true -use_fast_math $(GENCODE_FLAGS)
MPI_CFLAGS = -I${MPIROOT}/include
FFTW_CFLAGS = -I${FFTWROOT}/include

ifeq ($(OS),linux)
CUDA_LFLAGS = -L$(CUDAROOT)/lib64
else
ifeq ($(OS),titan)
CUDA_LFLAGS = -L$(CUDAROOT)/lib64
else
CUDA_LFLAGS = -L$(CUDAROOT)/lib
endif
endif
CUDA_LFLAGS += -lcudart -lcufft -lnvToolsExt

FFTW_LFLAGS = -L$(FFTWROOT)/lib -lfftw3 -lfftw3f

ifeq ($(MPI_FOUND), $(YES))
BINARIES = cpu_gse
endif
ifeq ($(CUDA_COMPILER), $(YES))
BINARIES += gpu_conv
endif

all: $(BINARIES)

cpu_gse : $(OBJS_GSE)
	$(CL) $(OPENMP_OPT) $(FFTW_LFLAGS) -o cpu_gse $(OBJS_GSE)

gpu_conv : $(OBJS_GPU_CONV)
	$(CL) $(OPENMP_OPT) $(CUDA_LFLAGS) -o gpu_conv $(OBJS_GPU_CONV)

clean: 
	rm -f *.o
	rm -f *.d
	rm -f *~
	rm -f $(BINARIES)

# Pull in dependencies that already exist
-include $(OBJS:.o=.d)

%.o : %.cu
	nvcc -c $(MPI_CFLAGS) $(NVCC_CFLAGS) $(DEFS) $<
	nvcc -M $(MPI_CFLAGS) $(NVCC_CFLAGS) $(DEFS) $*.cu > $*.d

%.o : %.cpp
	$(CC) -c $(CUDA_CFLAGS) $(FFTW_CFLAGS) $(DEFS) $<
	$(CC) -MM $(CUDA_CFLAGS) $(FFTW_CFLAGS) $(DEFS) $*.cpp > $*.d
