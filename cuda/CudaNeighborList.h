#ifndef CUDANEIGHBORLIST_H
#define CUDANEIGHBORLIST_H
//
// Neighbor list class
//
// (c) Antti-Pekka Hynninen 2014
// aphynninen@hotmail.com
//
#include <iostream>
#include <vector>
#include <cuda.h>
#include "CudaTopExcl.h"
#include "CudaNeighborListSort.h"
#include "CudaNeighborListBuild.h"

template <int tilesize>
class CudaNeighborList {

private:
  // Number of zones
  int numZone;

  // Number of lists
  int numList;

  // Zone parameters
  int h_ZoneParam_len;
  ZoneParam_t* h_ZoneParam;

  int d_ZoneParam_len;
  ZoneParam_t* d_ZoneParam;

  // Topological exclusions
  const CudaTopExcl& topExcl;

  // List sorting
  std::vector< CudaNeighborListSort* > sorter;

  // List building
  std::vector< CudaNeighborListBuild<tilesize>* > builder;

  // Total number of columns for each list
  //std::vector<int> ncol_tot;

  // NlistParam lists
  std::vector<NlistParam_t*> d_NlistParam;
  std::vector<NlistParam_t*> h_NlistParam;

  // Events
  std::vector<cudaEvent_t> build_event;
  cudaEvent_t glo2loc_reset_event;

  // Image boundaries: -1, 0, 1
  int imx_lo, imx_hi;
  int imy_lo, imy_hi;
  int imz_lo, imz_hi;

  // Number of z-cells in each column
  int col_ncellz_len;
  int *col_ncellz;

  // Starting cell index for each column
  int col_cell_len;
  int* col_cell;

  // Index sorting performed by this class
  // Usage: ind_sorted[i] = j : Atom j belongs to position i
  int ind_sorted_len;
  int *ind_sorted;

  // Atom indices where each cell start
  int cell_patom_len;
  int *cell_patom;

  // (icellx, icelly, icellz, izone) for each cell
  int cell_xyz_zone_len;
  int4 *cell_xyz_zone;

  // Cell z-boundaries
  int cell_bz_len;
  float *cell_bz;

  // Approximate upper bound for number of cells
  int ncell_max;

  // Bounding boxes
  int bb_len;
  bb_t *bb;

  // Flag for testing neighborlist build
  bool test;

  void set_NlistParam(cudaStream_t stream);
  void get_NlistParam();

  void sort_realloc(const int indList);

public:
  CudaNeighborList(const CudaTopExcl& topExcl,
		   const int nx, const int ny, const int nz);
  ~CudaNeighborList();

  void registerList(std::vector<int>& numIntZone,
		    std::vector< std::vector<int> >& intZones,
		    const char* filename=NULL);

  void sort(const int indList, const int *zone_patom,
	    const float4 *xyzq,
	    float4 *xyzq_sorted,
	    int *loc2glo,
	    cudaStream_t stream=0);

  void build(const int indList, const int* zone_patom,
	     const float boxx, const float boxy, const float boxz,
	     const float rcut,
	     const float4 *xyzq, const int *loc2glo,
	     cudaStream_t stream=0);

  //void reset();
  
  int *get_glo2loc() {return topExcl.get_glo2loc();}

  int *get_ind_sorted() {return ind_sorted;}

  int getNumList() {return numList;}

  void set_test(bool test_in) {
    this->test = test_in;
    for (int indList=0;indList < numList;indList++) {
      sorter.at(indList)->set_test(test_in);
      builder.at(indList)->set_test(test_in);
    }
  }

  CudaNeighborListBuild<tilesize>& getBuilder(const int indList) {
    assert(indList >= 0 && indList < numList);
    return *builder.at(indList);
  }

  void analyze() {
    for (int indList=0;indList < numList;indList++) {
      builder.at(indList)->analyze();
    }
  }

};

#endif // CUDANEIGHBORLIST_H
