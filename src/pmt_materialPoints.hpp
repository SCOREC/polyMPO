#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "pmt_utils.hpp"
#include "pmt_mesh.hpp"

#include <pumipic_kktypes.hpp>
#include <particle_structs.hpp>
#include <Kokkos_Core.hpp>


namespace polyMpmTest{

using particle_structs::DPS;
using particle_structs::MemberTypes;

typedef bool mp_status_t;
typedef int  mp_elem_id_t;
typedef double mp_vec2d_t[2];
typedef double mp_vec3d_t[3];
typedef double mp_sym_mat3d_t[6];
typedef bool mp_flag_basis_t;
typedef double mp_basis_t[maxVtxsPerElm];
typedef double mp_basis_grad2d_t[maxVtxsPerElm][2];
typedef double mp_constv_mdl_param_t[12];

enum MaterialPointSlice {
  MP_STATUS = 0,
  MP_CUR_ELM_ID = 1,
  MP_TGT_ELM_ID = 2,
  MP_CUR_POSITION_LAT_LON = 3,
  MP_TGT_POSITION_LAT_LON = 4,
  MP_CUR_POSITION_XYZ = 5,
  MP_TGT_POSITION_XYZ = 6,
  MP_FLAG_BASIS_VALS = 7,
  MP_BASIS_VALS = 8,
  MP_BASIS_GRAD_VALS = 9,
  MP_VELOCITY = 10,
  MP_STRAIN_RATE = 11,
  MP_STRESS = 12,
  MP_STRESS_DIV = 13,
  SHEAR_TRACTION = 14,
  MP_CONST_MDL_PARAM = 15
};

typedef MemberTypes<mp_elem_id_t,mp_elem_id_t,mp_elem_id_t,mp_vec2d_t,mp_vec2d_t,mp_vec3d_t,mp_vec3d_t,mp_elem_id_t,mp_basis_t,mp_basis_grad2d_t,mp_vec2d_t,mp_sym_mat3d_t,mp_sym_mat3d_t,mp_vec3d_t,mp_vec3d_t,mp_constv_mdl_param_t> MaterialPointTypes;
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
    auto getPositions() { return getData<MP_CUR_POSITION_XYZ>(); }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

