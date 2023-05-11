#include "pmt_materialPoints.hpp"

namespace polyMpmTest {

PS* createDPS(int numElms, int numMPs, Vector2View positions, IntView mpsPerElm, IntView mp2elm) {
  PS::kkGidView elmGids("elementGlobalIds", numElms); //TODO - initialize this to [0..numElms)
  auto mpInfo = ps::createMemberViews<MaterialPointTypes>(numMPs);
  auto mpPositions = ps::getMemberView<MaterialPointTypes, 0>(mpInfo);
  Kokkos::parallel_for("setMPinfo", numMPs, KOKKOS_LAMBDA(int i) {
    mpPositions(i,0) = positions[i][0];
    mpPositions(i,1) = positions[i][1];
  });
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,Kokkos::AUTO);
  return new DPS<MaterialPointTypes>(policy, numElms, numMPs, mpsPerElm, elmGids, mp2elm, mpInfo);
}

}
