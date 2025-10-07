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

        auto dualTriangleAreaEntry = meshFields2TypeAndString.at(MeshF_DualTriangleArea);
        PMT_ALWAYS_ASSERT(dualTriangleAreaEntry.first == MeshFType_VtxBased);
    }
    
    void Mesh::setMeshElmBasedFieldSize(){
        PMT_ALWAYS_ASSERT(meshEdit_);
        
        auto elmMassMapEntry = meshFields2TypeAndString.at(MeshF_ElmMass);
        PMT_ALWAYS_ASSERT(elmMassMapEntry.first == MeshFType_ElmBased);
        elmMass_      = MeshFView<MeshF_ElmMass>(elmMassMapEntry.second,numElms_);
                
        elmCenterXYZ_ = MeshFView<MeshF_ElmCenterXYZ>(meshFields2TypeAndString.at(MeshF_ElmCenterXYZ).second, numElms_);

        elmCenterGnomProj_= MeshFView<MeshF_ElmCenterGnomProj>(meshFields2TypeAndString.at(MeshF_ElmCenterGnomProj).second, numElms_);

        vtxGnomProj_ = MeshFView<MeshF_VtxGnomProj>(meshFields2TypeAndString.at(MeshF_VtxGnomProj).second, numElms_);
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

    KOKKOS_INLINE_FUNCTION
    void computeGnomonicProjectionAtPoint( const Vec3d& vtxCoord, const double gnomProjElmCenter[4], double& outX, double& outY) {
      const double iDen = 1.0 / ( gnomProjElmCenter[1] * gnomProjElmCenter[3] * vtxCoord[0] +
                                  gnomProjElmCenter[0] * gnomProjElmCenter[3] * vtxCoord[1] +
                                  gnomProjElmCenter[2] * vtxCoord[2]);
      outX = iDen * (vtxCoord[1] * gnomProjElmCenter[1] -
                    vtxCoord[1] * gnomProjElmCenter[0]);
      outY = iDen * (vtxCoord[2] * gnomProjElmCenter[3] - vtxCoord[1] * gnomProjElmCenter[2] * gnomProjElmCenter[0] -
                     vtxCoord[0] * gnomProjElmCenter[1] * gnomProjElmCenter[2]);
    }

    void Mesh::setGnomonicProjection(bool isRotated, double radius){
      std::cout<<__FUNCTION__<<std::endl;
      auto gnomProjVtx = getMeshField<MeshF_VtxGnomProj>();
      auto gnomProjElmCenter = getMeshField<MeshF_ElmCenterGnomProj>();
      
      auto vtxCoords  = getMeshField<MeshF_VtxCoords>();
      auto elmCenters = getMeshField<MeshF_ElmCenterXYZ>();
      auto elm2VtxConn = getElm2VtxConn();  


      Kokkos::parallel_for("setGnomprojCenter", numElms_, KOKKOS_LAMBDA(const int iElm){
        
        Vec3d elmCenter(elmCenters(iElm, 0), elmCenters(iElm, 1), elmCenters(iElm, 2));
        if(isRotated){
        
        }
        auto cos2LatR = elmCenter[0]*elmCenter[0] + elmCenter[1]*elmCenter[1];
        auto invR = 1.0/ sqrt(cos2LatR + elmCenter[2]*elmCenter[2]);
        auto cosLatR = sqrt(cos2LatR);

        gnomProjElmCenter(iElm, 0) = elmCenter[1]/cosLatR;
        gnomProjElmCenter(iElm, 1) = elmCenter[0]/cosLatR; 
        gnomProjElmCenter(iElm, 2) = invR*elmCenter[2];
        gnomProjElmCenter(iElm, 3) = invR*cosLatR;

        int nVtxE = elm2VtxConn(iElm,0);
        for(int i=0; i<nVtxE; i++){
          int vID = elm2VtxConn(iElm, i+1) - 1;
          Vec3d vtxCord(vtxCoords(vID, 0), vtxCoords(vID, 1), vtxCoords(vID, 2));
          if (isRotated){

          }

          double outX, outY;
          computeGnomonicProjectionAtPoint(vtxCord, &gnomProjElmCenter(iElm, 0), outX, outY);

          gnomProjVtx(iElm, i, 0) = outX;
          gnomProjVtx(iElm, i, 1) = outY;
        }
      });
    }
   
} // namespace polyMPO
