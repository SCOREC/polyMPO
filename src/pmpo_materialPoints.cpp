#include "pmpo_materialPoints.hpp"

namespace polyMPO {

namespace {
template<typename MemSpace = defaultSpace, typename View>
pumipic::MemberTypeViews createInternalMemberViews(int numMPs, View mp2elm, View mpAppID){
  auto mpInfo = ps::createMemberViews<MaterialPointTypes, MemSpace>(numMPs);
  auto mpCurElmPos_m = ps::getMemberView<MaterialPointTypes, MPF_Cur_Elm_ID, MemSpace>(mpInfo);
  auto mpAppID_m = ps::getMemberView<MaterialPointTypes, MPF_MP_APP_ID, MemSpace>(mpInfo);
  auto mpStatus_m = ps::getMemberView<MaterialPointTypes, MPF_Status, MemSpace>(mpInfo);
  auto policy = Kokkos::RangePolicy<typename MemSpace::execution_space>(typename MemSpace::execution_space(), 0, numMPs);
  Kokkos::parallel_for("setMPinfo", policy, KOKKOS_LAMBDA(int i) {
    mpCurElmPos_m(i) = mp2elm(i);
    mpStatus_m(i) = MP_ACTIVE;
    mpAppID_m(i) = mpAppID(i);
  });
  return mpInfo;
}

PS* createDPS(int numElms, int numMPs, DoubleVec3dView positions, IntView mpsPerElm, IntView mp2elm) {
  PS::kkGidView elmGids("elementGlobalIds", numElms); //TODO - initialize this to [0..numElms)
  auto mpInfo = ps::createMemberViews<MaterialPointTypes>(numMPs);
  auto mpPositions = ps::getMemberView<MaterialPointTypes, MPF_Cur_Pos_XYZ>(mpInfo);
  auto mpCurElmPos = ps::getMemberView<MaterialPointTypes, MPF_Cur_Elm_ID>(mpInfo);
  auto mpStatus = ps::getMemberView<MaterialPointTypes, MPF_Status>(mpInfo);
  auto mpAppID_m = ps::getMemberView<MaterialPointTypes, MPF_MP_APP_ID>(mpInfo);
  Kokkos::parallel_for("setMPinfo", numMPs, KOKKOS_LAMBDA(int i) {
    mpPositions(i,0) = positions(i,0);
    mpPositions(i,1) = positions(i,1);
    mpPositions(i,2) = positions(i,2);
    mpCurElmPos(i) = mp2elm(i);
    mpStatus(i) = MP_ACTIVE; 
    mpAppID_m(i) = i;
  });
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,Kokkos::AUTO);
  auto dps = new DPS<MaterialPointTypes>(policy, numElms, numMPs, mpsPerElm, elmGids, mp2elm, mpInfo);
  ps::destroyViews<MaterialPointTypes>(mpInfo);
  return dps;
}

PS* createDPS(int numElms, int numMPs, IntView mpsPerElm, IntView mp2elm, IntView mpAppID) {
  PS::kkGidView elmGids("elementGlobalIds", numElms); //TODO - initialize this to [0..numElms)
  auto mpInfo = createInternalMemberViews(numMPs, mp2elm, mpAppID);
  Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,Kokkos::AUTO);
  auto dps = new DPS<MaterialPointTypes>(policy, numElms, numMPs, mpsPerElm, elmGids, mp2elm, mpInfo);
  ps::destroyViews<MaterialPointTypes>(mpInfo);
  return dps;
}
}

MaterialPoints::MaterialPoints(int numElms, int numMPs, DoubleVec3dView positions, IntView mpsPerElm, IntView mp2elm) {
  MPs = createDPS(numElms, numMPs, positions, mpsPerElm, mp2elm);
  maxAppID = numMPs; //this ctor does not support inactive MPs
  operating_mode = MP_RELEASE;
};

MaterialPoints::MaterialPoints(int numElms, int numMPs, IntView mpsPerElm, IntView mp2elm, IntView mpAppID) {
  MPs = createDPS(numElms, numMPs, mpsPerElm, mp2elm, mpAppID);
  updateMaxAppID();
  operating_mode = MP_RELEASE;
};

MaterialPoints::~MaterialPoints() {
  if(MPs != nullptr)
    delete MPs;
}

void MaterialPoints::startRebuild(IntView tgtElm, int newNumMPs, IntView newMP2elm, IntView newMPAppID) {
  rebuildNumNewMPs = newNumMPs;
  rebuildNewMP2elm = newMP2elm;
  rebuildtgtElm = tgtElm;
  auto newMP2elm_h = Kokkos::create_mirror_view_and_copy(hostSpace(), newMP2elm);
  auto newMPAppID_h = Kokkos::create_mirror_view_and_copy(hostSpace(), newMPAppID);
  buildSlices = createInternalMemberViews<hostSpace>(newNumMPs, newMP2elm_h, newMPAppID_h);
}

void MaterialPoints::rebuild(IntView tgtElm, int newNumMPs, IntView newMP2elm, IntView newMPAppID) {
  startRebuild(tgtElm, newNumMPs, newMP2elm, newMPAppID);
  finishRebuild();
}
}