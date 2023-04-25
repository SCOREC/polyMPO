#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "pmt_utils.hpp"

#include <pumipic_kktypes.hpp>
#include <particle_structs.hpp>
#include <Kokkos_Core.hpp>


namespace polyMpmTest{

using particle_structs::DPS;
using particle_structs::MemberTypes;

typedef MemberTypes<double[2], int, bool> MaterialPointTypes;
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
      //const int numElms = 1; //FIXME - putting all MPs in one element
      PS::kkLidView ppe("particlesPerElement", 1); //FIXME - putting all MPs in one element
      Kokkos::deep_copy(Kokkos::subview(ppe,0),count);
      PS::kkGidView elmGids("elementGlobalIds", count);
      //Kokkos::TeamPolicy<ExeSpace> policy(num_elems,32);
      //particleStructure = new ps::DPS<Types,Memspace>(policy, numElms, count, ppe,
      //  element_gids, particle_elements, particle_info);
    };

    int getCount() { return count_; }
    Vector2View getPositions() { return positions_; }
    BoolView isActive() { return isActive_; }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

