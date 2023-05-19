#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "pmt_utils.hpp"

#include <pumipic_kktypes.hpp>
#include <particle_structs.hpp>
#include <Kokkos_Core.hpp>


namespace polyMpmTest{

using particle_structs::DPS;
using particle_structs::MemberTypes;

typedef std::int32_t LO; //local ordinal - count things on a single process that are less than 2^31
typedef std::int64_t GO; //global ordinal - count things across multiple processes that are greater than 2^31
typedef double Real2[2];
typedef double Real3[3];

enum MaterialPointSlice {
  MP_POSITION_R3 = 0
};

typedef MemberTypes<Real2> MaterialPointTypes;
typedef ps::ParticleStructure<MaterialPointTypes> PS;

PS* createDPS(int numElms, int numMPs, Vector2View positions, IntView mpsPerElm, IntView mp2elm);

class MaterialPoints {
  private:
    PS* MPs;

  public:
    MaterialPoints() : MPs(nullptr) {};
    MaterialPoints(int numElms, int numMPs, Vector2View positions, IntView mpsPerElm, IntView mp2elm) {
      MPs = createDPS(numElms, numMPs, positions, mpsPerElm, mp2elm);
    };
    ~MaterialPoints() {
      if(MPs != nullptr)
        delete MPs;
    }
    void rebuild(IntView materialPoints2Elm) {
      assert(materialPoints2Elm.size() == static_cast<size_t>(MPs->nPtcls()));
      if( materialPoints2Elm.size() < static_cast<size_t>(MPs->capacity()) ) {
        Kokkos::resize(Kokkos::WithoutInitializing, materialPoints2Elm, MPs->capacity());
      }
      MPs->rebuild(materialPoints2Elm);
    }
    template <int index>
    auto getData() {
      return MPs->get<index>();
    }
    template <typename FunctorType>
    void parallel_for(FunctorType kernel, std::string name="") {
      ps::parallel_for(MPs, kernel, name);
    }
    int getCount() { return MPs->nPtcls(); }
    auto getPositions() { return getData<0>(); }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

