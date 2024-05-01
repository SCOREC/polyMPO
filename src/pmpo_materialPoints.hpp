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
using hostSpace = Kokkos::HostSpace;
using defaultSpace = Kokkos::DefaultExecutionSpace::memory_space;

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
typedef std::function<int()> IntFunc;

enum MaterialPointSlice {
  MPF_Status = 0,
  MPF_Cur_Elm_ID,
  MPF_Tgt_Elm_ID,
  MPF_Cur_Pos_Rot_Lat_Lon,
  MPF_Tgt_Pos_Rot_Lat_Lon,
  MPF_Cur_Pos_XYZ,
  MPF_Tgt_Pos_XYZ,
  MPF_Flag_Basis_Vals,
  MPF_Basis_Vals,
  MPF_Basis_Grad_Vals,
  MPF_Mass,
  MPF_Vel,
  MPF_Rot_Lat_Lon_Incr,
  MPF_Strain_Rate,
  MPF_Stress,
  MPF_Stress_Div,
  MPF_Shear_Traction,
  MPF_Constv_Mdl_Param,
  MPF_MP_APP_ID
};

enum Operating_Mode{
  MP_RELEASE,
  MP_DEBUG
};

template <MaterialPointSlice> const int mpSliceToMeshFieldSize;
template <> const int mpSliceToMeshFieldSize < MPF_Status              > = 1;
template <> const int mpSliceToMeshFieldSize < MPF_Cur_Elm_ID          > = 0;
template <> const int mpSliceToMeshFieldSize < MPF_Tgt_Elm_ID          > = 0;
template <> const int mpSliceToMeshFieldSize < MPF_Cur_Pos_Rot_Lat_Lon > = 2;
template <> const int mpSliceToMeshFieldSize < MPF_Tgt_Pos_Rot_Lat_Lon > = 2;
template <> const int mpSliceToMeshFieldSize < MPF_Cur_Pos_XYZ         > = 3;
template <> const int mpSliceToMeshFieldSize < MPF_Tgt_Pos_XYZ         > = 3;
template <> const int mpSliceToMeshFieldSize < MPF_Flag_Basis_Vals     > = 1;
template <> const int mpSliceToMeshFieldSize < MPF_Basis_Vals          > = maxVtxsPerElm;
template <> const int mpSliceToMeshFieldSize < MPF_Basis_Grad_Vals     > = maxVtxsPerElm*2;
template <> const int mpSliceToMeshFieldSize < MPF_Mass                > = 1;
template <> const int mpSliceToMeshFieldSize < MPF_Vel                 > = 2;
template <> const int mpSliceToMeshFieldSize < MPF_Rot_Lat_Lon_Incr    > = 2;
template <> const int mpSliceToMeshFieldSize < MPF_Strain_Rate         > = 6;
template <> const int mpSliceToMeshFieldSize < MPF_Stress              > = 6;
template <> const int mpSliceToMeshFieldSize < MPF_Stress_Div          > = 3;
template <> const int mpSliceToMeshFieldSize < MPF_Shear_Traction      > = 3;
template <> const int mpSliceToMeshFieldSize < MPF_Constv_Mdl_Param    > = 12;
template <> const int mpSliceToMeshFieldSize < MPF_MP_APP_ID           > = 1;

const static std::vector<std::pair<MaterialPointSlice, MaterialPointSlice>>
        mpSliceSwap = {{MPF_Cur_Elm_ID, MPF_Tgt_Elm_ID},
                       {MPF_Cur_Pos_Rot_Lat_Lon, MPF_Tgt_Pos_Rot_Lat_Lon},
                       {MPF_Cur_Pos_XYZ, MPF_Tgt_Pos_XYZ}};

typedef MemberTypes<mp_flag_t,              //MPF_Status
                    mp_elm_id_t,            //MPF_Cur_Elm_ID
                    mp_elm_id_t,            //MPF_Tgt_Elm_ID
                    mp_vec2d_t,             //MPF_Cur_Pos_Rot_Lat_Lon
                    mp_vec2d_t,             //MPF_Tgt_Pos_Rot_Lat_Lon
                    mp_vec3d_t,             //MPF_Cur_Pos_XYZ
                    mp_vec3d_t,             //MPF_Tgt_Pos_XYZ
                    mp_flag_t,              //MPF_Flag_Basis_Vals
                    mp_basis_t,             //MPF_Basis_Vals
                    mp_basis_grad2d_t,      //MPF_Basis_Grad_Vals
                    mp_sclr_t,              //MPF_Mass //TODO: test Mass in assembly
                    mp_vec2d_t,             //MPF_Vel
                    mp_vec2d_t,             //MPF_Rot_Lat_Lon_Incr
                    mp_sym_mat3d_t,         //MPF_Strain_Rate
                    mp_sym_mat3d_t,         //MPF_Stress
                    mp_vec3d_t,             //MPF_Stress_Div
                    mp_vec3d_t,             //MPF_Shear_Traction
                    mp_constv_mdl_param_t,  //MPF_Constv_Mdl_Param
                    mp_id_t                 //MPF_APP_ID
                    >MaterialPointTypes;
typedef ps::ParticleStructure<MaterialPointTypes> PS;

struct RebuildHelper {
  bool ongoing = false;
  int addedNumMPs;
  pumipic::MemberTypeViews addedSlices_h;
  IntView addedMP2elm;
  IntView allTgtElm;
  Kokkos::View<const int*, hostSpace> addedMPMask_h;
};

class MaterialPoints {
  private:
    PS* MPs;
    int elmIDoffset = -1;
    int maxAppID = -1;
    bool isRotatedFlag = false;
    Operating_Mode operating_mode;
    RebuildHelper rebuildFields;
    IntFunc getAppID;

  public:
    MaterialPoints() : MPs(nullptr) {};
    MaterialPoints(int numElms, int numMPs, DoubleVec3dView positions, IntView mpsPerElm, IntView mp2elm);
    MaterialPoints(int numElms, int numMPs, IntView mpsPerElm, IntView mp2elm, IntView mpAppID);
    ~MaterialPoints();

    void rebuild(IntView addedMP2elm, IntView addedMPAppID);
    void startRebuild(IntView tgtElm, int addedNumMPs, IntView addedMP2elm, IntView addedMPAppID, Kokkos::View<const int*> addedMPMask);
    void finishRebuild();
    bool rebuildOngoing();
    
    template<int mpSliceIndex, typename mpSliceData>
    typename std::enable_if<mpSliceData::rank==1>::type
    setRebuildMPSlice(mpSliceData mpSliceIn);

    template<int mpSliceIndex, typename mpSliceData>
    typename std::enable_if<mpSliceData::rank==2>::type
    setRebuildMPSlice(mpSliceData mpSliceIn);

    void setAppIDFunc(IntFunc getAppIDIn);
    int getNextAppID();

    void rebuild() {
      IntView tgtElm("tgtElm", MPs->capacity());
      auto tgtMpElm = MPs->get<MPF_Tgt_Elm_ID>();
      auto setTgtElm = PS_LAMBDA(const int&, const int& mp, const int& mask) {
        if(mask) {
          tgtElm(mp) = tgtMpElm(mp);
        }
      };
      ps::parallel_for(MPs, setTgtElm, "setTargetElement");
      MPs->rebuild(tgtElm);
    }
    void updateMPElmID(){
      auto curElmID = MPs->get<MPF_Cur_Elm_ID>();
      auto tgtElmID = MPs->get<MPF_Tgt_Elm_ID>();
      auto swap = PS_LAMBDA(const int&, const int& mp, const int& mask) {
        if(mask){
            curElmID(mp) = tgtElmID(mp);
            tgtElmID(mp) = -1;
        }
      };
      ps::parallel_for(MPs, swap, "swap");
    }
    void updateMaxAppID() {
      auto mpAppID_m = MPs->get<MPF_MP_APP_ID>();
      maxAppID = 0;
      Kokkos::parallel_reduce("setMax" , MPs->nPtcls(),
        KOKKOS_LAMBDA(const int i, int & valueToUpdate) {
          if ( mpAppID_m(i) > valueToUpdate ) valueToUpdate = mpAppID_m(i);
        },
        Kokkos::Max<int>(maxAppID)
      );
    }
    template <MaterialPointSlice mpfIndexCur, MaterialPointSlice mpfIndexTgt>
    void updateMPSlice(){
      auto curData = MPs->get<mpfIndexCur>();
      auto tgtData = MPs->get<mpfIndexTgt>();
      const int numEntriesCur = mpSliceToMeshFieldSize<mpfIndexCur>;
      const int numEntriesTgt = mpSliceToMeshFieldSize<mpfIndexTgt>;
      PMT_ALWAYS_ASSERT(numEntriesCur == numEntriesTgt);
      
      auto swap = PS_LAMBDA(const int&, const int& mp, const int& mask) {
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
        updateMPSlice<MPF_Cur_Pos_Rot_Lat_Lon,MPF_Tgt_Pos_Rot_Lat_Lon>();
        updateMPSlice<MPF_Cur_Pos_XYZ,MPF_Tgt_Pos_XYZ>();
    }
    void updateRotLatLonAndXYZ2Tgt(const double radius){
        auto curPosRotLatLon = MPs->get<MPF_Cur_Pos_Rot_Lat_Lon>();
        auto tgtPosRotLatLon = MPs->get<MPF_Tgt_Pos_Rot_Lat_Lon>();
        auto tgtPosXYZ = MPs->get<MPF_Tgt_Pos_XYZ>();
        auto rotLatLonIncr = MPs->get<MPF_Rot_Lat_Lon_Incr>();
        
        auto updateRotLatLon = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
            if(mask){
                auto rotLat = curPosRotLatLon(mp,0) + rotLatLonIncr(mp,0); // phi
                auto rotLon = curPosRotLatLon(mp,1) + rotLatLonIncr(mp,1); // lambda
                auto geoLat = rotLat;
                auto geoLon = rotLon;
                tgtPosRotLatLon(mp,0) = geoLat;
                tgtPosRotLatLon(mp,1) = geoLon;
                // x = cosLon cosLat, y = sinLon cosLat, z = sinLat
                tgtPosXYZ(mp,0) = radius * std::cos(geoLon) * std::cos(geoLat);
                tgtPosXYZ(mp,1) = radius * std::sin(geoLon) * std::cos(geoLat);
                tgtPosXYZ(mp,2) = radius * std::sin(geoLat); 
            } 
        };
        if(isRotatedFlag){
            //TODO rotation lat lon calc
            fprintf(stderr, "rotational lat lon in MP is not support yet!");
            PMT_ALWAYS_ASSERT(false);
        } 
        ps::parallel_for(MPs, updateRotLatLon,"updateRotationalLatitudeLongitude"); 
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
    int getNumElems() { return MPs->nElems(); }
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
    bool getRotatedFlag() {
      return isRotatedFlag;
    }
    void setRotatedFlag(bool flagSet) {
      isRotatedFlag = flagSet;
    }

    // MUTATOR  
    template <MaterialPointSlice index> void fillData(double value);//use PS_LAMBDA fill up to 1
};// End MaterialPoints

template <MaterialPointSlice index>
void MaterialPoints::fillData(double value){
  auto mpData = getData<index>();
  const int numEntries = mpSliceToMeshFieldSize<index>;
  auto setValue = PS_LAMBDA(const int&, const int& mp, const int& mask){
    if(mask) { //if material point is 'active'/'enabled'
      for(int i=0; i<numEntries; i++){
          mpData(mp,i) = value;
      }
    }
  };
  parallel_for(setValue, "setValue");
}

template<int mpSliceIndex, typename mpSliceData>
typename std::enable_if<mpSliceData::rank==1>::type
MaterialPoints::setRebuildMPSlice(mpSliceData mpSliceIn) {
  auto mpSliceIn_h = Kokkos::create_mirror_view_and_copy(hostSpace(), mpSliceIn);
  auto mpSliceAdded_h = ps::getMemberView<MaterialPointTypes, mpSliceIndex, hostSpace>(rebuildFields.addedSlices_h);
  auto mpAppIDAdded_h = ps::getMemberView<MaterialPointTypes, MPF_MP_APP_ID, hostSpace>(rebuildFields.addedSlices_h);
  for (int i=0; i < mpSliceAdded_h.extent(0); i++)
    mpSliceAdded_h(i) = mpSliceIn_h(mpAppIDAdded_h(i));
}

template<int mpSliceIndex, typename mpSliceData>
typename std::enable_if<mpSliceData::rank==2>::type
MaterialPoints::setRebuildMPSlice(mpSliceData mpSliceIn) {
  auto mpSliceIn_h = Kokkos::create_mirror_view_and_copy(hostSpace(), mpSliceIn);
  auto mpSliceAdded_h = ps::getMemberView<MaterialPointTypes, mpSliceIndex, hostSpace>(rebuildFields.addedSlices_h);
  auto mpAppIDAdded_h = ps::getMemberView<MaterialPointTypes, MPF_MP_APP_ID, hostSpace>(rebuildFields.addedSlices_h);
  for (int i=0; i < mpSliceAdded_h.extent(0); i++)
  for (int j=0; j < mpSliceAdded_h.extent(1); j++)
    mpSliceAdded_h(i,j) = mpSliceIn_h(j,mpAppIDAdded_h(i));
}

}
#endif

