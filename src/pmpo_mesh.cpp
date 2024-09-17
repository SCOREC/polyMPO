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

    void Mesh::setMeshVtxBasedFieldSize(){
        PMT_ALWAYS_ASSERT(meshEdit_);

        auto vtxCoordsMapEntry = meshFields2TypeAndString.at(MeshF_VtxCoords);
        PMT_ALWAYS_ASSERT(vtxCoordsMapEntry.first == MeshFType_VtxBased);
        vtxCoords_ = MeshFView<MeshF_VtxCoords>(vtxCoordsMapEntry.second,numVtxs_);

        auto vtxRotLatMapEntry = meshFields2TypeAndString.at(MeshF_VtxRotLat);
        PMT_ALWAYS_ASSERT(vtxRotLatMapEntry.first == MeshFType_VtxBased);
        vtxRotLat_ = MeshFView<MeshF_VtxRotLat>(vtxRotLatMapEntry.second,numVtxs_);

        auto vtxVelMapEntry = meshFields2TypeAndString.at(MeshF_Vel);
        PMT_ALWAYS_ASSERT(vtxVelMapEntry.first == MeshFType_VtxBased);
        vtxVel_ = MeshFView<MeshF_Vel>(vtxVelMapEntry.second,numVtxs_);

        auto vtxMassMapEntry = meshFields2TypeAndString.at(MeshF_VtxMass);
        PMT_ALWAYS_ASSERT(vtxMassMapEntry.first == MeshFType_VtxBased);
        vtxMass_ = MeshFView<MeshF_VtxMass>(vtxMassMapEntry.second,numVtxs_);

        auto vtxOnSurfVeloIncrMapEntry = meshFields2TypeAndString.at(MeshF_OnSurfVeloIncr);
        PMT_ALWAYS_ASSERT(vtxOnSurfVeloIncrMapEntry.first == MeshFType_VtxBased);
        vtxOnSurfVeloIncr_ = MeshFView<MeshF_OnSurfVeloIncr>(vtxOnSurfVeloIncrMapEntry.second,numVtxs_);

        auto vtxOnSurfDispIncrMapEntry = meshFields2TypeAndString.at(MeshF_OnSurfDispIncr);
        PMT_ALWAYS_ASSERT(vtxOnSurfDispIncrMapEntry.first == MeshFType_VtxBased);
        vtxOnSurfDispIncr_ = MeshFView<MeshF_OnSurfDispIncr>(vtxOnSurfDispIncrMapEntry.second,numVtxs_);

        auto vtxRotLatLonIncrMapEntry = meshFields2TypeAndString.at(MeshF_RotLatLonIncr);
        PMT_ALWAYS_ASSERT(vtxRotLatLonIncrMapEntry.first == MeshFType_VtxBased);
        vtxRotLatLonIncr_ = MeshFView<MeshF_RotLatLonIncr>(vtxRotLatLonIncrMapEntry.second,numVtxs_);
    }
    
    void Mesh::setMeshElmBasedFieldSize(){
        PMT_ALWAYS_ASSERT(meshEdit_);
        
        auto elmMassMapEntry = meshFields2TypeAndString.at(MeshF_ElmMass);
        PMT_ALWAYS_ASSERT(elmMassMapEntry.first == MeshFType_ElmBased);
        elmMass_      = MeshFView<MeshF_ElmMass>(elmMassMapEntry.second,numElms_);
        elmCenterXYZ_ = MeshFView<MeshF_ElmCenterXYZ>(meshFields2TypeAndString.at(MeshF_ElmCenterXYZ).second,numElms_);
    }

    void Mesh::computeRotLatLonIncr(){
        Kokkos::Timer timer;
        PMT_ALWAYS_ASSERT(geomType_ == geom_spherical_surf);
       
        auto dispIncr = getMeshField<MeshF_OnSurfDispIncr>();
        auto rotLatLonIncr = getMeshField<MeshF_RotLatLonIncr>();
        auto lat = getMeshField<MeshF_VtxRotLat>();
        auto sphereRadius = getSphereRadius();
        PMT_ALWAYS_ASSERT(sphereRadius > 0); 
        Kokkos::parallel_for("set nEdgesPerElm", numVtxs_, KOKKOS_LAMBDA(const int iVtx){
            // Lat [iVtx,0] = dispIncrY [iVtx,1] /R
            // Lon [iVtx,1] = dispIncrX [iVtx,0] /(R*cos(lat))
            rotLatLonIncr(iVtx, 0) = dispIncr(iVtx, 1)/sphereRadius;
            rotLatLonIncr(iVtx, 1) = dispIncr(iVtx, 0)/(sphereRadius * std::cos(lat(iVtx)));
        });
        pumipic::RecordTime("PolyMPO_computeRotLatLonIncr", timer.seconds());
    }

} // namespace polyMPO
