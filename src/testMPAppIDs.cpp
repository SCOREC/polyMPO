#ifndef TESTMPAPPIDS_H
#define TESTMPAPPIDS_H

#include "pmpo_createTestMPMesh.hpp"
#include "pmpo_c.h"

using space_t = Kokkos::DefaultExecutionSpace::memory_space;

template<typename DataT>
using kkViewHostU = Kokkos::View<
          DataT,
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

typedef kkViewHostU<double*> kkDblViewHostU;//TODO:put it to mesh.hpp             
typedef kkViewHostU<double*[vec2d_nEntries]> kkVec2dViewHostU;//TODO:put it to mesh.hpp
typedef kkViewHostU<double**> kkDbl2dViewHostU;//TODO:put it somewhere else (maybe)
typedef kkViewHostU<int**> kkInt2dViewHostU;//TODO:put it somewhere else (maybe)
typedef kkViewHostU<int*> kkIntViewHostU;//TODO:put it somewhere else (maybe)

template <typename DataT>
auto create_mirror_view_and_copy(DataT array, const int size){
  kkViewHostU<DataT> temp_host(array, size);
  return Kokkos::create_mirror_view_and_copy(space_t(), temp_host);
}

extern "C" {
  void polympo_testFortranPointer_f(MPMesh_ptr p_mpmesh, const int numMPs, const int* allMP2Elm);
}

void polympo_testFortranPointer_f(MPMesh_ptr p_mpmesh, 
                              const int numMPs, // total number of MPs which is GREATER than or equal to number of active MPs
                              const int* allMP2Elm) {
//   checkMPMeshValid(p_mpmesh);
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getCount());
  PMT_ALWAYS_ASSERT(numMPs >= p_MPs->getMaxAppID());

  int internalMPCapacity = p_MPs->getCapacity(); // pumipic expects full capacity to rebuild
  auto mp2Elm = create_mirror_view_and_copy(allMP2Elm, internalMPCapacity);
  Kokkos::View<int*> mp2Elm_d("mp2Elm_d", internalMPCapacity);
  Kokkos::deep_copy(mp2Elm_d, mp2Elm);

  int numAddedMPs = 2;
  Kokkos::View<int*> added_mp2Elm_d("added_mp2Elm_d", numAddedMPs);
  Kokkos::parallel_for("set addedMP2Elm", numAddedMPs, KOKKOS_LAMBDA (const int i) {
    added_mp2Elm_d(i) = i;
  });
  
  Kokkos::View<int*> addedMPMask_d("addedMPMask_d", numMPs + numAddedMPs);
  Kokkos::parallel_for("set addedMPMask", numMPs + numAddedMPs, KOKKOS_LAMBDA (const int i) {
    if (i >= numMPs) addedMPMask_d(i) = 1;
    else addedMPMask_d(i) = 0;
  });

  std::vector<int> added_mpIDs(numAddedMPs);
  for(int i=0; i<numAddedMPs; i++)
    added_mpIDs[i] = p_MPs->getNextAppID();
  auto added_mpIDs_d = create_mirror_view_and_copy(added_mpIDs.data(), numAddedMPs);

  p_MPs->startRebuild(mp2Elm_d, numAddedMPs, added_mp2Elm_d, added_mpIDs_d, addedMPMask_d);
  p_MPs->finishRebuild();

  auto newAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::parallel_for("print APP ID", numAddedMPs, KOKKOS_LAMBDA (const int i) {
    printf(" Added %d\n", newAppID(numMPs+i));
    assert(added_mpIDs_d[i] == newAppID(numMPs+i));
  });
}

#endif