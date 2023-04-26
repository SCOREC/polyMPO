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
    int numElms_;
    int numMPs_;
    Vector2View positions_;
    BoolView isActive_;
    PS* MPs;

  public:
    MaterialPoints() : MPs(nullptr) {};
    MaterialPoints(int numElms, int numMPs, Vector2View positions, BoolView isActive):
                    numElms_(numElms), numMPs_(numMPs), positions_(positions), isActive_(isActive){
      PS::kkLidView ppe("particlesPerElement", numElms_);
      PS::kkGidView elmGids("elementGlobalIds", numElms_);
      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms_,32);
      MPs = new DPS<MaterialPointTypes>(policy, numElms_, numMPs_, ppe, elmGids);
    };
    ~MaterialPoints() {
      if(MPs != nullptr)
        delete MPs;
    }

    int getCount() { return numMPs_; }
    Vector2View getPositions() { return positions_; }
    BoolView isActive() { return isActive_; }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

