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
} // namespace polyMPO
