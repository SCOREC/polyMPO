#include "pmpo_mesh.hpp"

namespace polyMPO{
    bool Mesh::checkMeshType(int meshType){
        for(const auto& validType : validMeshType){
            if(meshType == validType){
               return true;
            }
        }
        return false; 
    }

    bool Mesh::checkGeomType(int geomType){
        for(const auto& validType : validGeomType){
            if(geomType == validType){
                return true;
            }
        }
        return false;
    }
} // namespace polyMPO
