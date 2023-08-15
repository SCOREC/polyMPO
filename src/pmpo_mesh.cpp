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

    void Mesh::setMeshFieldSize(int numVtxs){
        PMT_ALWAYS_ASSERT(meshEdit_);
        vtxVel_ = DoubleVec2dView(meshFields2String.at(MeshF_Vel).second,numVtxs);
        vtxOnSurfVeloIncr_ = DoubleVec2dView(meshFields2String.at(MeshF_OnSurfVeloIncr).second,numVtxs);
        vtxOnSurfDispIncr_ = DoubleVec2dView(meshFields2String.at(MeshF_OnSurfDispIncr).second,numVtxs);
    }

} // namespace polyMPO
