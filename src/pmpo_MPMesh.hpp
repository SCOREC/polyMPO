#ifndef POLYMPO_MPM_H
#define POLYMPO_MPM_H

#include "pmpo_utils.hpp"
#include "pmpo_mesh.hpp"
#include "pmpo_materialPoints.hpp"

namespace polyMPO{

template <MaterialPointSlice> const MeshFieldIndex mpSliceToMeshFieldIndex;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Status              > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Cur_Elm_ID          > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Tgt_Elm_ID          > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Cur_Pos_Rot_Lat_Lon > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Tgt_Pos_Rot_Lat_Lon > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Cur_Pos_XYZ         > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Tgt_Pos_XYZ         > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Flag_Basis_Vals     > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Basis_Vals          > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Basis_Grad_Vals     > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Mass                > = MeshF_Unsupported;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Vel                 > = MeshF_Vel;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Rot_Lat_Lon_Incr    > = MeshF_RotLatLonIncr;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Strain_Rate         > = MeshF_Unsupported;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Stress              > = MeshF_Unsupported;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Stress_Div          > = MeshF_Unsupported;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Shear_Traction      > = MeshF_Unsupported;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_Constv_Mdl_Param    > = MeshF_Invalid;
template <> const MeshFieldIndex mpSliceToMeshFieldIndex < MPF_MP_APP_ID           > = MeshF_Invalid;

#define maxMPsPerElm 8

class MPMesh{
  public:
    Mesh* p_mesh;
    MaterialPoints* p_MPs;

    int reconsOption;
    std::map<MaterialPointSlice, std::function<void()>> reconstructSlice = std::map<MaterialPointSlice, std::function<void()>>();
    
    MPMesh(Mesh* inMesh, MaterialPoints* inMPs):
        p_mesh(inMesh), p_MPs(inMPs) {
    };
    ~MPMesh() {
      delete p_mesh;
      delete p_MPs;
    }

    void CVTTrackingEdgeCenterBased(Vec2dView dx);
    void CVTTrackingElmCenterBased(const int printVTPIndex = -1);
    void T2LTracking(Vec2dView dx);
    void push();

    DoubleView assemblyV0();
    template <MaterialPointSlice index>
    DoubleView wtScaAssembly();
    template <MaterialPointSlice index>
    Vec2dView wtVec2Assembly();
    template <MaterialPointSlice mpfIndex>
    void assembly(bool basisWeightFlag, bool massWeightFlag);

    template<MaterialPointSlice mpSliceIndex>
    void setReconstructSlice();
    void reconstructSlices();

    void printVTP_mesh(int printVTPIndex);
};

}//namespace polyMPO end

#endif

