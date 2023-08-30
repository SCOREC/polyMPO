#include "pmpo_materialPoints.hpp"

namespace polyMPO {

PS* createDPS(int numElms, int numMPs, DoubleVec3dView positions, IntView mpsPerElm, IntView mp2elm) {
  PS::kkGidView elmGids("elementGlobalIds", numElms); //TODO - initialize this to [0..numElms)
  auto mpInfo = ps::createMemberViews<MaterialPointTypes>(numMPs);
  auto mpPositions = ps::getMemberView<MaterialPointTypes, MPF_Cur_Pos_XYZ>(mpInfo);
  auto mpCurElmPos = ps::getMemberView<MaterialPointTypes, MPF_Cur_Elm_ID>(mpInfo);
  auto mpStatus = ps::getMemberView<MaterialPointTypes, MPF_Status>(mpInfo);
  Kokkos::parallel_for("setMPinfo", numMPs, KOKKOS_LAMBDA(int i) {
    mpPositions(i,0) = positions(i,0);
    mpPositions(i,1) = positions(i,1);
    mpPositions(i,2) = positions(i,2);
    mpCurElmPos(i) = mp2elm(i);
    mpStatus(i) = 1; 
  });
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,Kokkos::AUTO);
  auto dps = new DPS<MaterialPointTypes>(policy, numElms, numMPs, mpsPerElm, elmGids, mp2elm, mpInfo);
  ps::destroyViews<MaterialPointTypes>(mpInfo);
  return dps;
}

PS* createDPS(int numElms, int numMPs, IntView mpsPerElm, IntView mp2elm) {
  PS::kkGidView elmGids("elementGlobalIds", numElms); //TODO - initialize this to [0..numElms)
  auto mpInfo = ps::createMemberViews<MaterialPointTypes>(numMPs);
  auto mpCurElmPos = ps::getMemberView<MaterialPointTypes, MPF_Cur_Elm_ID>(mpInfo);
  auto mpStatus = ps::getMemberView<MaterialPointTypes, MPF_Status>(mpInfo);
  Kokkos::parallel_for("setMPinfo", numMPs, KOKKOS_LAMBDA(int i) {
    mpCurElmPos(i) = mp2elm(i);
    mpStatus(i) = 1;
  });
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,Kokkos::AUTO);
  auto dps = new DPS<MaterialPointTypes>(policy, numElms, numMPs, mpsPerElm, elmGids, mp2elm, mpInfo);
  ps::destroyViews<MaterialPointTypes>(mpInfo);
  return dps;
}

}
