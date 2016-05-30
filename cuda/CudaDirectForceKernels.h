#ifndef CUDADIRECTFORCEKERNELS_H
#define CUDADIRECTFORCEKERNELS_H

#include <cuda.h>
#include "Bonded_struct.h"
#include "CudaNeighborListBuild.h"
#include "CudaBlock.h"
#include "CudaDirectForceTypes.h"
#include "CudaEnergyVirial.h"

#ifndef USE_TEXTURE_OBJECTS
texture<float2, 1, cudaReadModeElementType>* get_vdwparam_texref();
texture<float2, 1, cudaReadModeElementType>* get_vdwparam14_texref();
texture<float, 1, cudaReadModeElementType>* getBlockParamTexRef();
bool get_vdwparam_texref_bound();
bool get_vdwparam14_texref_bound();
bool getBlockParamTexRefBound();
void set_vdwparam_texref_bound(const bool val);
void set_vdwparam14_texref_bound(const bool val);
void setBlockParamTexRefBound(const bool val);
#endif

const int tilesize = 32;


template <typename AT, typename CT>
void calcForceKernelChoice(const int nblock_tot_in, const int nthread, const int shmem_size, cudaStream_t stream,
			   const int vdw_model, const int elec_model, const bool calc_energy, const bool calc_virial,
			   const CudaNeighborListBuild<32>& nlist,
			   const float* vdwparam, const int nvdwparam, const int* vdwtype,
#ifdef USE_TEXTURE_OBJECTS
			   cudaTextureObject_t& vdwParamTexObj,
#endif
			   const float4* xyzq, const int stride, AT* force,
			   Virial_t *virial, double *energy_vdw, double *energy_elec,
			   CudaBlock* cudaBlock=NULL, AT* biflam=NULL, AT* biflam2=NULL);

template <typename AT, typename CT>
void calcForce14KernelChoice(const int nblock, const int nthread, const int shmem_size, cudaStream_t stream,
			     const int vdw_model, const int elec_model, const bool calc_energy, const bool calc_virial,
			     const int nin14list, const xx14list_t* in14list, const int nex14list, const xx14list_t* ex14list,
			     const int nin14block, const int* vdwtype, const float* vdwparam14,
#ifdef USE_TEXTURE_OBJECTS
			     cudaTextureObject_t& vdwParam14TexObj,
#endif
			     const float4* xyzq, const float fscale, const int stride, AT* force,
			     Virial_t *virial, double *energy_vdw, double *energy_elec, double *energy_excl);

void updateDirectForceSetup(const DirectSettings_t* h_setup);

#endif //CUDADIRECTFORCEKERNELS_H
