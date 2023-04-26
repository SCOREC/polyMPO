#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "pmt_utils.hpp"

#include <pumipic_kktypes.hpp>
#include <particle_structs.hpp>
#include <Kokkos_Core.hpp>


namespace polyMpmTest{

using particle_structs::DPS;
using particle_structs::MemberTypes;
using pumipic::Vector3d;

typedef MemberTypes<Vector3d, int, int> MaterialPointTypes;
typedef ps::ParticleStructure<MaterialPointTypes> PS;

class MaterialPoints {
  private:
    Vector2View positions_;
    BoolView isActive_;
    PS* MPs;

  public:
    MaterialPoints() : MPs(nullptr) {};
    MaterialPoints(int numElms, int numMPs, Vector2View positions, BoolView isActive):
                    positions_(positions), isActive_(isActive){
      PS::kkLidView ppe("particlesPerElement", numElms);
      PS::kkGidView elmGids("elementGlobalIds", numElms);
      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,32);
      MPs = new DPS<MaterialPointTypes>(policy, numElms, numMPs, ppe, elmGids);
    };
    ~MaterialPoints() {
      if(MPs != nullptr)
        delete MPs;
    }
    void rebuild(IntView materialPoints2Elm) {
      assert(materialPoints2Elm.size() == MPs->nPtcls());
      if( materialPoints2Elm.size() < static_cast<size_t>(MPs->capacity()) ) {
        Kokkos::resize(Kokkos::WithoutInitializing, materialPoints2Elm, MPs->capacity());
      }
      MPs->rebuild(materialPoints2Elm);
    }
    int getCount() { return MPs->nPtcls(); }
    Vector2View getPositions() { return positions_; }
    BoolView isActive() { return isActive_; }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

