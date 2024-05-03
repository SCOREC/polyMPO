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

template <typename Type, int Size> struct MPTypeHelper { using type = Type[Size]; static const int size = Size; };
template <typename Type> struct MPTypeHelper <Type, 0> { using type = int; static const int size = 0; };

template <MaterialPointSlice> struct mpSliceToMeshField;
template <> struct mpSliceToMeshField < MPF_Status              > { using typeHelper = MPTypeHelper<int,     0>; };
template <> struct mpSliceToMeshField < MPF_Cur_Elm_ID          > { using typeHelper = MPTypeHelper<int,     0>; };
template <> struct mpSliceToMeshField < MPF_Tgt_Elm_ID          > { using typeHelper = MPTypeHelper<int,     0>; };
template <> struct mpSliceToMeshField < MPF_Cur_Pos_Rot_Lat_Lon > { using typeHelper = MPTypeHelper<double,  2>; };
template <> struct mpSliceToMeshField < MPF_Tgt_Pos_Rot_Lat_Lon > { using typeHelper = MPTypeHelper<double,  2>; };
template <> struct mpSliceToMeshField < MPF_Cur_Pos_XYZ         > { using typeHelper = MPTypeHelper<double,  3>; };
template <> struct mpSliceToMeshField < MPF_Tgt_Pos_XYZ         > { using typeHelper = MPTypeHelper<double,  3>; };
template <> struct mpSliceToMeshField < MPF_Flag_Basis_Vals     > { using typeHelper = MPTypeHelper<int,     0>; };
template <> struct mpSliceToMeshField < MPF_Basis_Vals          > { using typeHelper = MPTypeHelper<double,  maxVtxsPerElm>; };
template <> struct mpSliceToMeshField < MPF_Basis_Grad_Vals     > { using typeHelper = MPTypeHelper<double,  maxVtxsPerElm*2>; };
template <> struct mpSliceToMeshField < MPF_Mass                > { using typeHelper = MPTypeHelper<double,  1>; };
template <> struct mpSliceToMeshField < MPF_Vel                 > { using typeHelper = MPTypeHelper<double,  2>; };
template <> struct mpSliceToMeshField < MPF_Rot_Lat_Lon_Incr    > { using typeHelper = MPTypeHelper<double,  2>; };
template <> struct mpSliceToMeshField < MPF_Strain_Rate         > { using typeHelper = MPTypeHelper<double,  6>; };
template <> struct mpSliceToMeshField < MPF_Stress              > { using typeHelper = MPTypeHelper<double,  6>; };
template <> struct mpSliceToMeshField < MPF_Stress_Div          > { using typeHelper = MPTypeHelper<double,  3>; };
template <> struct mpSliceToMeshField < MPF_Shear_Traction      > { using typeHelper = MPTypeHelper<double,  3>; };
template <> struct mpSliceToMeshField < MPF_Constv_Mdl_Param    > { using typeHelper = MPTypeHelper<double,  12>; };
template <> struct mpSliceToMeshField < MPF_MP_APP_ID           > { using typeHelper = MPTypeHelper<int,     0>; };

template <MaterialPointSlice slice> using mpSliceToMeshFieldType = typename mpSliceToMeshField<slice>::typeHelper::type;
template <MaterialPointSlice slice> const int mpSliceToMeshFieldSize = mpSliceToMeshField<slice>::typeHelper::size;

const static std::vector<std::pair<MaterialPointSlice, MaterialPointSlice>>
        mpSliceSwap = {{MPF_Cur_Elm_ID, MPF_Tgt_Elm_ID},
                       {MPF_Cur_Pos_Rot_Lat_Lon, MPF_Tgt_Pos_Rot_Lat_Lon},
                       {MPF_Cur_Pos_XYZ, MPF_Tgt_Pos_XYZ}};

typedef MemberTypes<mpSliceToMeshFieldType < MPF_Status              >,
                    mpSliceToMeshFieldType < MPF_Cur_Elm_ID          >,
                    mpSliceToMeshFieldType < MPF_Tgt_Elm_ID          >,
                    mpSliceToMeshFieldType < MPF_Cur_Pos_Rot_Lat_Lon >,
                    mpSliceToMeshFieldType < MPF_Tgt_Pos_Rot_Lat_Lon >,
                    mpSliceToMeshFieldType < MPF_Cur_Pos_XYZ         >,
                    mpSliceToMeshFieldType < MPF_Tgt_Pos_XYZ         >,
                    mpSliceToMeshFieldType < MPF_Flag_Basis_Vals     >,
                    mpSliceToMeshFieldType < MPF_Basis_Vals          >,
                    mpSliceToMeshFieldType < MPF_Basis_Grad_Vals     >,
                    mpSliceToMeshFieldType < MPF_Mass                >,
                    mpSliceToMeshFieldType < MPF_Vel                 >,
                    mpSliceToMeshFieldType < MPF_Rot_Lat_Lon_Incr    >,
                    mpSliceToMeshFieldType < MPF_Strain_Rate         >,
                    mpSliceToMeshFieldType < MPF_Stress              >,
                    mpSliceToMeshFieldType < MPF_Stress_Div          >,
                    mpSliceToMeshFieldType < MPF_Shear_Traction      >,
                    mpSliceToMeshFieldType < MPF_Constv_Mdl_Param    >,
                    mpSliceToMeshFieldType < MPF_MP_APP_ID           >
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

