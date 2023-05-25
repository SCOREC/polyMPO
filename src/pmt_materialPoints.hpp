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

//typedef bool mp_flag_t;
typedef int mp_flag_t;
typedef int  mp_elem_id_t;
typedef double mp_vec2d_t[2];
typedef double mp_vec3d_t[3];
typedef double mp_sym_mat3d_t[6];
typedef double mp_basis_t[maxVtxsPerElm];
typedef double mp_basis_grad2d_t[maxVtxsPerElm][2];
typedef double mp_constv_mdl_param_t[12];

enum MaterialPointSlice {
  MP_STATUS = 0,
  MP_CUR_ELM_ID = 1,
  MP_TGT_ELM_ID = 2,
  MP_CUR_POS_LAT_LON = 3,
  MP_TGT_POS_LAT_LON = 4,
  MP_CUR_POS_XYZ = 5,
  MP_TGT_POS_XYZ = 6,
  MP_FLAG_BASIS_VALS = 7,
  MP_BASIS_VALS = 8,
  MP_BASIS_GRAD_VALS = 9,
  MP_VEL = 10,
  MP_STRAIN_RATE = 11,
  MP_STRESS = 12,
  MP_STRESS_DIV = 13,
  MP_SHEAR_TRACTION = 14,
  MP_CONSTV_MDL_PARAM = 15
};

typedef MemberTypes<mp_flag_t,              //MP_STATUS
                    mp_elem_id_t,           //MP_CUR_ELM_ID
                    mp_elem_id_t,           //MP_TGT_ELM_ID
                    mp_vec2d_t,             //MP_CUR_POS_LAT_LON
                    mp_vec2d_t,             //MP_TGT_POS_LAT_LON
                    mp_vec3d_t,             //MP_CUR_POS_XYZ
                    mp_vec3d_t,             //MP_TGT_POS_XYZ
                    mp_flag_t,              //MP_FLAG_BASIS_VALS
                    mp_basis_t,             //MP_BASIS_VALS
                    mp_basis_grad2d_t,      //MP_BASIS_GRAD_VALS
                    mp_vec2d_t,             //MP_VEL
                    mp_sym_mat3d_t,         //MP_STRAIN_RATE
                    mp_sym_mat3d_t,         //MP_STRESS
                    mp_vec3d_t,             //MP_STRESS_DIV
                    mp_vec3d_t,             //MP_SHEAR_TRACTION
                    mp_constv_mdl_param_t   //MP_CONSTV_MDL_PARAM
                    >MaterialPointTypes;
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
    auto getPositions() { return getData<MP_CUR_POS_XYZ>(); }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

