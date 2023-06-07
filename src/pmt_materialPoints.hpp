#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "pmt_utils.hpp"
#include "pmt_mesh.hpp"

#include <pumipic_kktypes.hpp>
#include <particle_structs.hpp>
#include <Kokkos_Core.hpp>
#include <map>

namespace polyMpmTest{

using particle_structs::DPS;
using particle_structs::MemberTypes;

//typedef bool mp_flag_t;
typedef int mp_flag_t;
typedef int  mp_elem_id_t;
//TODO: 
typedef double mp_vec2d_t[2];
typedef double mp_vec3d_t[3];
typedef double mp_sym_mat3d_t[6];
typedef double mp_basis_t[maxVtxsPerElm];
typedef double mp_basis_grad2d_t[maxVtxsPerElm][2];
typedef double mp_constv_mdl_param_t[12];

enum MaterialPointSlice {
  MP_Status = 0,
  MP_Cur_Elm_ID,        //1
  MP_Tgt_Elm_ID, 
  MP_Cur_Pos_Lat_Lon,
  MP_Tgt_Pos_Lat_Lon,
  MP_Cur_Pos_XYZ,       //5
  MP_Tgt_Pos_XYZ,       
  MP_Flag_Basis_Vals,
  MP_Basis_Vals,
  MP_Basis_Grad_Vals,
  //MP_Basis_Grad_Vals,
  MP_Mass,              //10
  MP_Vel,               
  MP_Strain_Rate,
  MP_Stress,
  MP_Stress_Div,
  MP_Shear_Traction,    //15
  MP_ConstV_MDL_Param
};

//const std::map<MaterialPointSlice,std::pair<int,MeshFields>> 
const std::map<MaterialPointSlice,std::pair<int,MeshFieldIndex>> 
      mpSlice2MeshField = {{MP_Status,          {1,meshFieldInvalid}},
                           {MP_Cur_Elm_ID,      {1,meshFieldInvalid}},
                           {MP_Tgt_Elm_ID,      {1,meshFieldInvalid}},
                           {MP_Cur_Pos_Lat_Lon, {2,meshFieldInvalid}},
                           {MP_Tgt_Pos_Lat_Lon, {2,meshFieldInvalid}},
                           {MP_Cur_Pos_XYZ,{3,meshFieldInvalid}},
                           {MP_Tgt_Pos_XYZ,{3,meshFieldInvalid}},
                           {MP_Flag_Basis_Vals,{1,meshFieldInvalid}},
                           {MP_Basis_Vals,{maxVtxsPerElm,meshFieldInvalid}},
                           {MP_Basis_Grad_Vals,{-1,meshFieldInvalid}},//XXX:2d array
                           {MP_Mass,{1,meshFieldUnsupported}},
                           {MP_Vel,{2,meshFieldVelocity}},
                           {MP_Strain_Rate,{6,meshFieldUnsupported}},
                           {MP_Stress,{6,meshFieldUnsupported}},
                           {MP_Stress_Div,{3,meshFieldUnsupported}},
                           {MP_Shear_Traction,{3,meshFieldUnsupported}},
                           {MP_ConstV_MDL_Param,{12,meshFieldUnsupported}}};

typedef MemberTypes<mp_flag_t,              //MP_Status
                    mp_elem_id_t,           //MP_Cur_ELM_ID
                    mp_elem_id_t,           //MP_Tgt_ELM_ID
                    mp_vec2d_t,             //MP_Cur_Pos_Lat_Lon
                    mp_vec2d_t,             //MP_Tgt_Pos_Lat_Lon
                    mp_vec3d_t,             //MP_Cur_Pos_XYZ
                    mp_vec3d_t,             //MP_Tgt_Pos_XYZ
                    mp_flag_t,              //MP_Flag_Basis_Vals
                    mp_basis_t,             //MP_Basis_Vals
                    mp_basis_grad2d_t,      //MP_Basis_Grad_Vals
                    double,                 //MP_Mass
                    mp_vec2d_t,             //MP_Vel
                    mp_sym_mat3d_t,         //MP_Strain_Rate
                    mp_sym_mat3d_t,         //MP_Stress
                    mp_vec3d_t,             //MP_Stress_Div
                    mp_vec3d_t,             //MP_Shear_Traction
                    mp_constv_mdl_param_t   //MP_ConstV_MDL_Param
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
    auto getPositions() { return getData<MP_Cur_Pos_XYZ>(); }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

