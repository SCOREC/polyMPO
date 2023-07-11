#include "pmo_c.h"
#include "pmo_createTestMPMesh.hpp"
#include <stdio.h>

namespace{
    std::vector<mpmesh> mpMeshes;////store the mpMeshes that is legal
}

void polympo_initialize() {
  printf("polympo_initialize c++\n");
  int isMpiInit;
  MPI_Initialized(&isMpiInit);
  PMT_ALWAYS_ASSERT(isMpiInit);
  Kokkos::initialize();
}

void polympo_finalize() {
  printf("polympo_finalize c++\n");
  Kokkos::finalize();
}

mpmesh polympo_createMpMesh() {
  auto mesh = polyMpmTest::initTestMesh(1);
  auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
  auto mps = polyMpmTest::initTestMPs(mesh, mpPerElement);
  auto mpMeshReturn = (mpmesh) new polyMpmTest::MPMesh(mesh,mps);
  mpMeshes.push_back(mpMeshReturn);
  return mpMeshReturn;
}

void polympo_deleteMpMesh(mpmesh mpMeshIn) {
  //chech mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());

  mpMeshes.erase(mpMeshInIter);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  delete mpMesh;
}

/**
 * Attention: this typedef is LayoutLeft, meaning that the first
 * index is the contiguous one. This matches the Fortran and GPU conventions for
 * allocations.
 */
typedef Kokkos::View<
          double*[vec2d_nEntries],
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkDblViewHostU;//TODO:put it to mesh.hpp

void polympo_setMPVelArray(mpmesh mpMeshIn, int size, double* array) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_setMPVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);
  auto MPs = mpMesh->MPs;
  auto mpVel = MPs->getData<polyMpmTest::MPF_Vel>();
 
  auto mpVelCopy = polyMpmTest::DoubleVec2DView("mpVelNewValue",size);
  //copy the host array to the device
  Kokkos::deep_copy(mpVelCopy,arrayHost);
  
  //modify the MP array with the mpVelCopy copied from the host array
  auto setVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask) { 
      for(int i=0; i<vec2d_nEntries; i++){
        mpVel(mp,i) = mpVelCopy(mp,i);
      }
    }
  };
  MPs->parallel_for(setVel, "setVel to array");
}

void polympo_setMeshVelArray(mpmesh mpMeshIn, int size, double* array) {
  //check mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_setMeshVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = mpMesh->getMesh();
  auto vtxField = mesh.getMeshField<polyMpmTest::MeshF_Vel>();

  //check the size
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*vec2d_nEntries)==vtxField.size());

  //copy the host array to the device
  Kokkos::deep_copy(vtxField,arrayHost);
}

void polympo_getMPVelArray(mpmesh mpMeshIn, int size, double* array) {
  //chech mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_getMPVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = mpMesh->getMesh();
  auto MPs = mpMesh->MPs;
  auto mpVel = MPs->getData<polyMpmTest::MPF_Vel>();
  auto mpVelCopy = polyMpmTest::DoubleVec2DView("copyOfMPVel",size);
    
  PMT_ALWAYS_ASSERT(static_cast<size_t>(MPs->getCount()*vec2d_nEntries)==mpVelCopy.size());

  //the pumipic 'slices' are not guaranteed to be contiguous, so we run a
  //pumipic parallel_for loop to fill a contigious array which will be copied to
  //a contiguous host array
  auto copyVel = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask) { 
      for(int i=0; i<vec2d_nEntries; i++){
        mpVelCopy(mp,i) = mpVel(mp,i);
      }
    }
  };
  MPs->parallel_for(copyVel, "copy mpVel to mpVelCopy");

  //copy the device array to the host
  Kokkos::deep_copy(arrayHost,mpVelCopy); 
}

void polympo_getMeshVelArray(mpmesh mpMeshIn, int size, double* array) {
  //chech mpMesh is valid
  auto mpMeshInIter = std::find(mpMeshes.begin(),mpMeshes.end(),mpMeshIn);
  PMT_ALWAYS_ASSERT(mpMeshInIter != mpMeshes.end());
  printf("polympo_getMeshVelArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = mpMesh->getMesh();
  auto vtxField = mesh.getMeshField<polyMpmTest::MeshF_Vel>();

  //check the size
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*vec2d_nEntries)==vtxField.size());

  //copy the device array to the host
  Kokkos::deep_copy(arrayHost, vtxField);
}

void polympo_setCommunicator(MPI_Fint fcomm){
    MPI_Comm comm = MPI_Comm_f2c(fcomm);
    int commSize;
    MPI_Comm_size(comm,&commSize);
    printf("polympo_setCommunicator with a communicator size: %d\n",commSize);
}
