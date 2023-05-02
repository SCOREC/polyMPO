#include "pmt_materialPoints.hpp"

namespace polyMpmTest {

PS* createDPS(int numElms, int numMPs, Vector2View positions, BoolView isActive) {
  PS::kkLidView ppe("particlesPerElement", numElms);
  PS::kkGidView elmGids("elementGlobalIds", numElms);
  PS::kkLidView mpElms("mpParentElement", numMPs);
  auto mpInfo = ps::createMemberViews<MaterialPointTypes>(numMPs);
  auto mpPositions = ps::getMemberView<MaterialPointTypes, 0>(mpInfo);
  Kokkos::parallel_for("setMPinfo", numMPs, KOKKOS_LAMBDA(int i) {
    mpPositions(i,0) = positions[i][0];
    mpPositions(i,1) = positions[i][1];
    //active particles are assigned to inactive particles are marked with elm id -1
    mpElms[i] = isActive[i] ? 0 : -1;
  });
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,32);
  return new DPS<MaterialPointTypes>(policy, numElms, numMPs, ppe, elmGids, mpElms, mpInfo);
}

}