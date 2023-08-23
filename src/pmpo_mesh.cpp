#include "pmpo_mesh.hpp"

namespace polyMPO{
    bool Mesh::checkMeshType(int meshType){
        return (meshType >mesh_unrecognized_lower &&
                meshType <mesh_unrecognized_upper); 
    }

    bool Mesh::checkGeomType(int geomType){
        return (geomType >mesh_unrecognized_lower &&
                geomType <mesh_unrecognized_upper);
    }

    void Mesh::setMeshVtxBasedFieldSize(int numVtxs){
        PMT_ALWAYS_ASSERT(meshEdit_);

        auto vtxCoordsMapEntry = meshFields2TypeAndString.at(MeshF_VtxCoords);
        PMT_ALWAYS_ASSERT(vtxCoordsMapEntry.first == MeshFType_VtxBased);
        vtxCoords_ = DoubleVec3dView(vtxCoordsMapEntry.second,numVtxs);

        auto vtxVelMapEntry = meshFields2TypeAndString.at(MeshF_Vel);
        PMT_ALWAYS_ASSERT(vtxVelMapEntry.first == MeshFType_VtxBased);
        vtxVel_ = DoubleVec2dView(vtxVelMapEntry.second,numVtxs);

        auto vtxOnSurfVeloIncrMapEntry = meshFields2TypeAndString.at(MeshF_OnSurfVeloIncr);
        PMT_ALWAYS_ASSERT(vtxOnSurfVeloIncrMapEntry.first == MeshFType_VtxBased);
        vtxOnSurfVeloIncr_ = DoubleVec2dView(vtxOnSurfVeloIncrMapEntry.second,numVtxs);

        auto vtxOnSurfDispIncrMapEntry = meshFields2TypeAndString.at(MeshF_OnSurfDispIncr);
        PMT_ALWAYS_ASSERT(vtxOnSurfDispIncrMapEntry.first == MeshFType_VtxBased);
        vtxOnSurfDispIncr_ = DoubleVec2dView(vtxOnSurfDispIncrMapEntry.second,numVtxs);
    }

} // namespace polyMPO
