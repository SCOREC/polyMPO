#ifndef POLYMPO_MATERIALPOINTS_H
#define POLYMPO_MATERIALPOINTS_H

#include "pmpo_utils.hpp"
#include "pmpo_mesh.hpp"
#include "pmpo_defines.h"

#include <pumipic_kktypes.hpp>
#include <particle_structs.hpp>
#include <Kokkos_Core.hpp>
#include <map>

namespace polyMPO{

using particle_structs::DPS;
using particle_structs::MemberTypes;

//typedef bool mp_flag_t;
typedef int mp_flag_t;
typedef int mp_id_t;
typedef int  mp_elm_id_t;
typedef double mp_sclr_t[1];//TODO
typedef vec2d_t mp_vec2d_t;
typedef double mp_vec3d_t[3];//TODO
typedef double mp_sym_mat3d_t[6];//TODO
typedef double mp_basis_t[maxVtxsPerElm];
typedef double mp_basis_grad2d_t[maxVtxsPerElm*2];
typedef double mp_constv_mdl_param_t[12];

enum MaterialPointSlice {
  MPF_Status = 0,
  MPF_Cur_Elm_ID,        //1
  MPF_Tgt_Elm_ID, 
  MPF_Cur_Pos_Lat_Lon,
  MPF_Tgt_Pos_Lat_Lon,
  MPF_Cur_Pos_XYZ,       //5
  MPF_Tgt_Pos_XYZ,       
  MPF_Flag_Basis_Vals,
  MPF_Basis_Vals,
  MPF_Basis_Grad_Vals,
  MPF_Mass,              //10
  MPF_Vel,               
  MPF_Strain_Rate,
  MPF_Stress,
  MPF_Stress_Div,
  MPF_Shear_Traction,    //15
  MPF_Constv_Mdl_Param,
  MPF_MP_APP_ID
};

enum Operating_Mode{
  MP_RELEASE,
  MP_DEBUG
};

const static std::map<MaterialPointSlice, std::pair<int,MeshFieldIndex>> 
      mpSlice2MeshFieldIndex = {{MPF_Status,     {1,MeshF_Invalid}},
                           {MPF_Cur_Elm_ID,      {0,MeshF_Invalid}},
                           {MPF_Tgt_Elm_ID,      {0,MeshF_Invalid}},
                           {MPF_Cur_Pos_Lat_Lon, {2,MeshF_Invalid}},
                           {MPF_Tgt_Pos_Lat_Lon, {2,MeshF_Invalid}},
                           {MPF_Cur_Pos_XYZ,     {3,MeshF_Invalid}},
                           {MPF_Tgt_Pos_XYZ,     {3,MeshF_Invalid}},
                           {MPF_Flag_Basis_Vals, {1,MeshF_Invalid}},
                           {MPF_Basis_Vals,      {maxVtxsPerElm,MeshF_Invalid}},
                           {MPF_Basis_Grad_Vals, {maxVtxsPerElm*2,MeshF_Invalid}},
                           {MPF_Mass,            {1,MeshF_Unsupported}},
                           {MPF_Vel,             {2,MeshF_Vel}},
                           {MPF_Strain_Rate,     {6,MeshF_Unsupported}},
                           {MPF_Stress,          {6,MeshF_Unsupported}},
                           {MPF_Stress_Div,      {3,MeshF_Unsupported}},
                           {MPF_Shear_Traction,  {3,MeshF_Unsupported}},
                           {MPF_Constv_Mdl_Param,{12,MeshF_Unsupported}},
                           {MPF_MP_APP_ID,       {1,MeshF_Invalid}}};

const static std::vector<std::pair<MaterialPointSlice, MaterialPointSlice>>
        mpSliceSwap = {{MPF_Cur_Elm_ID, MPF_Tgt_Elm_ID},
                       {MPF_Cur_Pos_Lat_Lon, MPF_Tgt_Pos_Lat_Lon},
                       {MPF_Cur_Pos_XYZ, MPF_Tgt_Pos_XYZ}};

typedef MemberTypes<mp_flag_t,              //MP_Status
                    mp_elm_id_t,            //MP_Cur_Elm_ID
                    mp_elm_id_t,            //MP_Tgt_Elm_ID
                    mp_vec2d_t,             //MP_Cur_Pos_Lat_Lon
                    mp_vec2d_t,             //MP_Tgt_Pos_Lat_Lon
                    mp_vec3d_t,             //MP_Cur_Pos_XYZ
                    mp_vec3d_t,             //MP_Tgt_Pos_XYZ
                    mp_flag_t,              //MP_Flag_Basis_Vals
                    mp_basis_t,             //MP_Basis_Vals
                    mp_basis_grad2d_t,      //MP_Basis_Grad_Vals
                    mp_sclr_t,              //MP_Mass //TODO: test Mass in assembly
                    mp_vec2d_t,             //MP_Vel
                    mp_sym_mat3d_t,         //MP_Strain_Rate
                    mp_sym_mat3d_t,         //MP_Stress
                    mp_vec3d_t,             //MP_Stress_Div
                    mp_vec3d_t,             //MP_Shear_Traction
                    mp_constv_mdl_param_t,  //MP_Constv_Mdl_Param
                    mp_id_t                 //MP_APP_ID
                    >MaterialPointTypes;
typedef ps::ParticleStructure<MaterialPointTypes> PS;


PS* createDPS(int numElms, int numMPs, DoubleVec3dView positions, IntView mpsPerElm, IntView mp2elm);
PS* createDPS(int numElms, int numMPs, IntView mpsPerElm, IntView mp2elm, IntView mpAppID);
void updateMaxAppID();

pumipic::MemberTypeViews createInternalMemberViews(int newNumMPs, IntView newMp2elm, IntView newMpAppID);

class MaterialPoints {
  private:
    PS* MPs;
    int elmIDoffset = -1;
    int maxAppID = -1;
    Operating_Mode operating_mode;

  public:
    MaterialPoints() : MPs(nullptr) {};
    MaterialPoints(int numElms, int numMPs, DoubleVec3dView positions, IntView mpsPerElm, IntView mp2elm) {
      MPs = createDPS(numElms, numMPs, positions, mpsPerElm, mp2elm);
      maxAppID = numMPs; //this ctor does not support inactive MPs
      operating_mode = MP_RELEASE;
    };
    MaterialPoints(int numElms, int numMPs, IntView mpsPerElm, IntView mp2elm, IntView mpAppID) {
      MPs = createDPS(numElms, numMPs, mpsPerElm, mp2elm, mpAppID);
      updateMaxAppID();
      operating_mode = MP_RELEASE;
    };
    ~MaterialPoints() {
      if(MPs != nullptr)
        delete MPs;
    }
    void rebuild();
    void rebuild(IntView tgtElm, int newNumMPs, IntView newMP2elm, IntView newMPAppID);
    void updateMPElmID();
    void updateMaxAppID();
    template <MaterialPointSlice mpfIndexCur, MaterialPointSlice mpfIndexTgt>
    void updateMPSlice(){
      auto curData = MPs->get<mpfIndexCur>();
      auto tgtData = MPs->get<mpfIndexTgt>();
      const int numEntriesCur = mpSlice2MeshFieldIndex.at(mpfIndexCur).first;
      const int numEntriesTgt = mpSlice2MeshFieldIndex.at(mpfIndexTgt).first;
      PMT_ALWAYS_ASSERT(numEntriesCur == numEntriesTgt);
      
      auto swap = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if(mask){
          for(int i=0; i<numEntriesCur; i++){
            curData(mp,i) = tgtData(mp,i);
            tgtData(mp,i) = -1;
          }
        }
      };
      ps::parallel_for(MPs, swap, "swap");
    }
    void updateMPSliceAll(){
        updateMPElmID();
        updateMPSlice<MPF_Cur_Pos_Lat_Lon,MPF_Tgt_Pos_Lat_Lon>();
        updateMPSlice<MPF_Cur_Pos_XYZ,MPF_Tgt_Pos_XYZ>();
    }

    template <int index>
    auto getData() {
      return MPs->get<index>();
    }
    template <typename FunctorType>
    void parallel_for(FunctorType kernel, std::string name="") {
      ps::parallel_for(MPs, kernel, name);
    }
    int getCapacity() { return MPs->capacity(); }
    int getCount() { return MPs->nPtcls(); }
    auto getPositions() { return getData<MPF_Cur_Pos_XYZ>(); }

    Operating_Mode getOpMode() { return operating_mode; }
    void setOpMode(Operating_Mode op_mode) {
      operating_mode = op_mode;
    }
    void setElmIDoffset(int offset) {
      PMT_ALWAYS_ASSERT(offset == 0 || offset == 1);
      elmIDoffset = offset;
    }
    int getElmIDoffset() {
      PMT_ALWAYS_ASSERT(elmIDoffset == 0 || elmIDoffset == 1);
      return elmIDoffset;
    }
    int getMaxAppID() {
      PMT_ALWAYS_ASSERT(maxAppID != -1);
      return maxAppID;
    }

//MUTATOR  
    template <MaterialPointSlice index> void fillData(double value);//use PS_LAMBDA fill up to 1
    void T2LTracking(Vec2dView dx);
    
};

template <MaterialPointSlice index>
void MaterialPoints::fillData(double value){
    auto mpData = getData<index>();
    const int numEntries = mpSlice2MeshFieldIndex.at(index).first;
    auto setValue = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
        if(mask) { //if material point is 'active'/'enabled'
            for(int i=0; i<numEntries; i++){
                mpData(mp,i) = value;
            }
        }
    };
    parallel_for(setValue, "setValue");
}


}
#endif

