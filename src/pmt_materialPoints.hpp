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
    int count_;
    Vector2View positions_;
    BoolView isActive_;
    PS* MPs;

  public:
    MaterialPoints(){};
    MaterialPoints( int count, Vector2View positions, BoolView isActive):
                    count_(count), positions_(positions), isActive_(isActive){
      const int numElms = 1; //FIXME - putting all MPs in one element
      PS::kkLidView ppe("particlesPerElement", numElms);
      PS::kkGidView elmGids("elementGlobalIds", numElms);
      Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace> policy(numElms,32);
      MPs = new DPS<MaterialPointTypes>(policy, numElms, count, ppe, elmGids);
    };
    ~MaterialPoints() {
      delete MPs;
    }

    int getCount() { return count_; }
    Vector2View getPositions() { return positions_; }
    BoolView isActive() { return isActive_; }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

