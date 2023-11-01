#include "pmpo_materialPoints.hpp"

namespace polyMPO {

pumipic::MemberTypeViews createInternalMemberViews(int numMPs, IntView mp2elm, IntView mpAppID){
  auto mpInfo = ps::createMemberViews<MaterialPointTypes>(numMPs);
  auto mpCurElmPos_m = ps::getMemberView<MaterialPointTypes, MPF_Cur_Elm_ID>(mpInfo);
  auto mpAppID_m = ps::getMemberView<MaterialPointTypes, MPF_MP_APP_ID>(mpInfo);
  auto mpStatus_m = ps::getMemberView<MaterialPointTypes, MPF_Status>(mpInfo);
  Kokkos::parallel_for("setMPinfo", numMPs, KOKKOS_LAMBDA(int i) {
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

void MaterialPoints::rebuild() {
  IntView tgtElm("tgtElm", MPs->capacity());
  auto tgtMpElm = MPs->get<MPF_Tgt_Elm_ID>();
  auto setTgtElm = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) {
      tgtElm(mp) = tgtMpElm(mp);
    }
  };
  ps::parallel_for(MPs, setTgtElm, "setTargetElement");
  MPs->rebuild(tgtElm);
}

void MaterialPoints::rebuild(IntView tgtElm, int newNumMPs, IntView newMP2elm, IntView newMPAppID) {
  auto newMPInfo = createInternalMemberViews(newNumMPs, newMP2elm, newMPAppID);
  MPs->rebuild(tgtElm, newMP2elm, newMPInfo);
  updateMaxAppID();
}

void MaterialPoints::updateMPElmID(){
  auto curElmID = MPs->get<MPF_Cur_Elm_ID>();
  auto tgtElmID = MPs->get<MPF_Tgt_Elm_ID>();
  auto swap = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask){
        curElmID(mp) = tgtElmID(mp);
        tgtElmID(mp) = -1;
    }
  };
  ps::parallel_for(MPs, swap, "swap");
}

void MaterialPoints::updateMaxAppID() {
  auto mpInfo = ps::createMemberViews<MaterialPointTypes>(MPs->nPtcls());
  auto mpAppID_m = ps::getMemberView<MaterialPointTypes, MPF_MP_APP_ID>(mpInfo);
  maxAppID = 0;
  Kokkos::parallel_reduce("setMax" , mpAppID_m.size(),
    KOKKOS_LAMBDA(const int i, int & valueToUpdate) {
      if ( mpAppID_m(i) > valueToUpdate ) valueToUpdate = mpAppID_m(i) ;
    },
    Kokkos::Max<int>(maxAppID)
  );
}
}