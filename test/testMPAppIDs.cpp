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
  assert(p_MPs->getNumElems() >= numAddedMPs);
  Kokkos::View<int*> added_mp2Elm_d("added_mp2Elm_d", numAddedMPs);
  Kokkos::parallel_for("set addedMP2Elm", numAddedMPs, KOKKOS_LAMBDA (const int i) {
    added_mp2Elm_d(i) = i;
  });

  std::vector<int> added_mpIDs(numAddedMPs);
  for(int i=0; i<numAddedMPs; i++)
    added_mpIDs[i] = p_MPs->getNextAppID();
  auto added_mpIDs_d = create_mirror_view_and_copy(added_mpIDs.data(), numAddedMPs);

  //Maps used to check values after rebuild
  Kokkos::UnorderedMap<int, int> added_MPs_data(numAddedMPs+1);
  Kokkos::parallel_for("set added_MPs_data", numAddedMPs, KOKKOS_LAMBDA (const int i) {
    added_MPs_data.insert(added_mpIDs_d(i), added_mp2Elm_d(i));
  });

  Kokkos::UnorderedMap<int, int> old_MPs_data(numAddedMPs+1);
  auto oldAppIDs = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  auto setOldAppIDs = PS_LAMBDA(const int& e, const int& mp, const int& mask) {
    if(mask)
      old_MPs_data.insert(oldAppIDs(mp), e);
  };
  p_MPs->parallel_for(setOldAppIDs, "setOldAppIDs");

  p_MPs->rebuild(added_mp2Elm_d, added_mpIDs_d);

  //Assert rebuild worked
  auto newAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();
  Kokkos::View<int*> numAddedMPsAfter("numAddedMPsAfter", 1);
  auto checkAddedMPs = PS_LAMBDA(const int& e, const int& mp, const int& mask) {
    if(mask) {
      if (added_MPs_data.exists(newAppID(mp))) {
        int index = added_MPs_data.find(newAppID(mp));
        Kokkos::atomic_increment(&numAddedMPsAfter(0));
        assert(e == added_MPs_data.value_at(index));
        added_MPs_data.insert(newAppID(mp), -1); //reset
      }
      else if (old_MPs_data.exists(newAppID(mp))) {
        int index = old_MPs_data.find(newAppID(mp));
        assert(e == old_MPs_data.value_at(index));
        old_MPs_data.insert(newAppID(mp), -1); //reset
      }
      else Kokkos::abort("Material point in wrong place!\n");
    }
  };
  p_MPs->parallel_for(checkAddedMPs, "checkAddedMPs");

  int numAddedMPsAfter_h = pumipic::getLastValue(numAddedMPsAfter);
  PMT_ALWAYS_ASSERT(numAddedMPsAfter_h == numAddedMPs);
}

#endif