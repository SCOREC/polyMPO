#ifndef TESTMPAPPIDS_H
#define TESTMPAPPIDS_H

#include "pmpo_createTestMPMesh.hpp"
#include "pmpo_defines.h"
#include "pmpo_c.h"

extern "C" {
  void testAppIDPointer(MPMesh_ptr p_mpmesh);
}

void testAppIDPointer(MPMesh_ptr p_mpmesh) {
  auto p_MPs = ((polyMPO::MPMesh*)p_mpmesh)->p_MPs;

  int numAddedMPs = 2;
  Kokkos::View<int*> added_mp2Elm_d("added_mp2Elm_d", numAddedMPs);
  Kokkos::parallel_for("set addedMP2Elm", numAddedMPs, KOKKOS_LAMBDA (const int i) {
    added_mp2Elm_d(i) = i;
  });

  std::vector<int> added_mpIDs(numAddedMPs);
  for(int i=0; i<numAddedMPs; i++)
    added_mpIDs[i] = p_MPs->getNextAppID();
  auto added_mpIDs_d = create_mirror_view_and_copy(added_mpIDs.data(), numAddedMPs);
  int prevNumMPs = p_MPs->getCount();

  p_MPs->rebuild(added_mp2Elm_d, added_mpIDs_d);

  //Assert rebuild worked
  auto newAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::parallel_for("print APP ID", numAddedMPs, KOKKOS_LAMBDA (const int i) {
    assert(added_mpIDs_d[i] == newAppID(prevNumMPs+i));
  });
}

#endif