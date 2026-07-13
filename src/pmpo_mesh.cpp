#include "pmpo_mesh.hpp"

namespace polyMPO{
  bool Mesh::checkMeshType(int meshType){
    return (meshType >mesh_unrecognized_lower && meshType <mesh_unrecognized_upper); 
  }

  bool Mesh::checkGeomType(int geomType){
    return (geomType >mesh_unrecognized_lower && geomType <mesh_unrecognized_upper);
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
    dualTriangleArea_ = MeshFView<MeshF_DualTriangleArea>(dualTriangleAreaEntry.second,numVtxs_);

    auto tanLatVertexRotOverRadiusEntry = meshFields2TypeAndString.at(MeshF_TanLatVertexRotatedOverRadius);
    PMT_ALWAYS_ASSERT(tanLatVertexRotOverRadiusEntry.first == MeshFType_VtxBased);
    tanLatVertexRotatedOverRadius_ = MeshFView<MeshF_TanLatVertexRotatedOverRadius>(tanLatVertexRotOverRadiusEntry.second, numVtxs_);

    auto interiorVertexEntry = meshFields2TypeAndString.at(MeshF_InteriorVertex);
    PMT_ALWAYS_ASSERT(interiorVertexEntry.first == MeshFType_VtxBased);
    interiorVertex_ = MeshFView<MeshF_InteriorVertex>(interiorVertexEntry.second, numVtxs_);

    uBdryVtxNormal_ = MeshFView<MeshF_uBdryVtxNormal>(meshFields2TypeAndString.at(MeshF_uBdryVtxNormal).second, numVtxs_, 2);

    vBdryVtxNormal_ = MeshFView<MeshF_vBdryVtxNormal>(meshFields2TypeAndString.at(MeshF_vBdryVtxNormal).second, numVtxs_, 2);

    auto stressDivergenceEntry = meshFields2TypeAndString.at(MeshF_StressDivergence);
    PMT_ALWAYS_ASSERT(stressDivergenceEntry.first == MeshFType_VtxBased);
    stressDivergence_ = MeshFView<MeshF_StressDivergence>(stressDivergenceEntry.second, numVtxs_, 2);

    auto solveVelocityEntry = meshFields2TypeAndString.at(MeshF_SolveVelocity);
    PMT_ALWAYS_ASSERT(solveVelocityEntry.first == MeshFType_VtxBased);
    solveVelocity_ = MeshFView<MeshF_SolveVelocity>(solveVelocityEntry.second, numVtxs_);

    auto totalMassVtxEntry = meshFields2TypeAndString.at(MeshF_TotalMassVtx);
    PMT_ALWAYS_ASSERT(totalMassVtxEntry.first == MeshFType_VtxBased);
    totalMassVtx_ = MeshFView<MeshF_TotalMassVtx>(totalMassVtxEntry.second, numVtxs_);

    //For grid solve, the source terms
    airStress_ = MeshFView<MeshF_AirStress>(meshFields2TypeAndString.at(MeshF_AirStress).second, numVtxs_);

    surfaceTilt_ = MeshFView<MeshF_SurfaceTilt>(meshFields2TypeAndString.at(MeshF_SurfaceTilt).second, numVtxs_);
    
    totalMassFVtx_ = MeshFView<MeshF_TotalMassFVtx>(meshFields2TypeAndString.at(MeshF_TotalMassFVtx).second, numVtxs_);
    
    oceanStress_ = MeshFView<MeshF_OceanStress>(meshFields2TypeAndString.at(MeshF_OceanStress).second, numVtxs_);
    
    oceanStressCoeff_ = MeshFView<MeshF_OceanStressCoeff>(meshFields2TypeAndString.at(MeshF_OceanStressCoeff).second, numVtxs_);

    oceanVelocity_ = MeshFView<MeshF_OceanVelocity>(meshFields2TypeAndString.at(MeshF_OceanVelocity).second, numVtxs_);

  }

  void Mesh::setMeshElmBasedFieldSize(){
    PMT_ALWAYS_ASSERT(meshEdit_);

    auto elmMassMapEntry = meshFields2TypeAndString.at(MeshF_ElmMass);
    PMT_ALWAYS_ASSERT(elmMassMapEntry.first == MeshFType_ElmBased);
    elmMass_      = MeshFView<MeshF_ElmMass>(elmMassMapEntry.second,numElms_);

    elmCenterXYZ_ = MeshFView<MeshF_ElmCenterXYZ>(meshFields2TypeAndString.at(MeshF_ElmCenterXYZ).second, numElms_);

    elmCenterGnomProj_= MeshFView<MeshF_ElmCenterGnomProj>(meshFields2TypeAndString.at(MeshF_ElmCenterGnomProj).second, numElms_);

    vtxGnomProj_ = MeshFView<MeshF_VtxGnomProj>(meshFields2TypeAndString.at(MeshF_VtxGnomProj).second, numElms_);

    solveStress_ = MeshFView<MeshF_SolveStress>(meshFields2TypeAndString.at(MeshF_SolveStress).second, numElms_);
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

  void Mesh::setGnomonicProjection(bool isRotated){
    std::cout<<__FUNCTION__<<std::endl;
    auto gnomProjVtx = getMeshField<MeshF_VtxGnomProj>();
    auto gnomProjElmCenter = getMeshField<MeshF_ElmCenterGnomProj>();

    auto vtxCoords  = getMeshField<MeshF_VtxCoords>();
    auto elmCenters = getMeshField<MeshF_ElmCenterXYZ>();
    auto elm2VtxConn = getElm2VtxConn();  

    Kokkos::parallel_for("setGnomprojCenter", numElms_, KOKKOS_LAMBDA(const int iElm){
      Vec3d elmCenter(elmCenters(iElm, 0), elmCenters(iElm, 1), elmCenters(iElm, 2));
      if(isRotated){
        elmCenter[0] = - elmCenters(iElm, 2);
        elmCenter[2] =   elmCenters(iElm, 0);
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
        if(isRotated){
          vtxCord[0] = - vtxCoords(vID, 2);
          vtxCord[2] = vtxCoords(vID, 0);
        }

        double outX, outY;
        auto gnomProjElmCenter_sub = Kokkos::subview(gnomProjElmCenter, iElm, Kokkos::ALL);
        computeGnomonicProjectionAtPoint(vtxCord, gnomProjElmCenter_sub, outX, outY);

        gnomProjVtx(iElm, i, 0) = outX;
        gnomProjVtx(iElm, i, 1) = outY;
      }
    });
  }

  void Mesh::calcOceanStressCoeff(){
    int numVerticesOwned = getNumVerticesOwned();
    auto iceAreaVtx = getMeshField<polyMPO::MeshF_VtxMass>();
    auto oceanStressCoeff = getMeshField<MeshF_OceanStressCoeff>();
    auto velocity = getMeshField<MeshF_Vel>();
    auto solve_velocity = getMeshField<MeshF_SolveVelocity>();
    auto oceanVelocity = getMeshField<MeshF_OceanVelocity>();

    Kokkos::parallel_for("calcOceanStressCoeff", numVerticesOwned, KOKKOS_LAMBDA(const int vtx){
      if(solve_velocity(vtx) == 0) return;
      auto relVelSq = pow(oceanVelocity(vtx, 0) - velocity(vtx, 0), 2) +  pow(oceanVelocity(vtx, 1) - velocity(vtx, 1), 2);
      oceanStressCoeff(vtx, 0) = 0.00536 * 1026.0 * iceAreaVtx(vtx, 0) * sqrt(relVelSq);
    });
  }

  void Mesh::gridSolveGPU(){
    //Mesh Fields
    int numVerticesOwned = getNumVerticesOwned(); 
    auto totalMassVtx = getMeshField<MeshF_TotalMassVtx>();
    auto totalMassFVtx = getMeshField<MeshF_TotalMassFVtx>();
    auto airStress = getMeshField<MeshF_AirStress>();
    auto surfaceTiltForce = getMeshField<MeshF_SurfaceTilt>();
    auto oceanStress = getMeshField<MeshF_OceanStress>();
    auto oceanStressCoeff = getMeshField<MeshF_OceanStressCoeff>();
    auto solve_velocity = getMeshField<MeshF_SolveVelocity>();
    auto stressDivergence = getMeshField<MeshF_StressDivergence>();
    auto velocity = getMeshField<MeshF_Vel>();
    auto elasticTimeStep = getElasticTimeStep(); 

    double sinOceanTurningAngle=0.0;
    double cosOceanTurningAngle=1.0;

    Kokkos::parallel_for("SolveGridVelocity", numVerticesOwned, KOKKOS_LAMBDA(const int vtx){
      if(solve_velocity(vtx) == 0) return;
      double a, b, c, d, s, rhs_u, rhs_v, denom ;
      
      s = (totalMassFVtx(vtx, 0) >= 0) ? 1.0 : -1.0;
      a = totalMassVtx(vtx, 0)/elasticTimeStep + oceanStressCoeff(vtx, 0) * cosOceanTurningAngle;
      b = -totalMassFVtx(vtx, 0) - oceanStressCoeff(vtx, 0) * sinOceanTurningAngle * s * totalMassFVtx(vtx, 0);
      c = totalMassFVtx(vtx, 0)  + oceanStressCoeff(vtx, 0) * sinOceanTurningAngle * s * totalMassFVtx(vtx, 0);
      d = totalMassVtx(vtx, 0)/elasticTimeStep + oceanStressCoeff(vtx, 0) * cosOceanTurningAngle; 
      rhs_u = stressDivergence(vtx, 0) + airStress(vtx, 0) + surfaceTiltForce(vtx, 0) + oceanStressCoeff(vtx, 0)*oceanStress(vtx, 0) +
              (totalMassVtx(vtx, 0) * velocity(vtx, 0))/elasticTimeStep;
      rhs_v = stressDivergence(vtx, 1) + airStress(vtx, 1) + surfaceTiltForce(vtx, 1) + oceanStressCoeff(vtx, 0)*oceanStress(vtx, 1) +
              (totalMassVtx(vtx, 0) * velocity(vtx, 1))/elasticTimeStep;
      denom = a*d - b*c;
      
      velocity(vtx, 0) = (d*rhs_u-b*rhs_v)/denom;
      velocity(vtx, 1) = (a*rhs_v-c*rhs_u)/denom;
    });
  }


  void Mesh::aggregateDeluDyn(){
    int numVtx = getNumVertices();
    auto elasticTimeStep = getElasticTimeStep();

    constexpr MeshFieldIndex mfIndex_surfDispIncr = MeshF_OnSurfDispIncr; 
    auto meshField_surfDispIncr = getMeshField<mfIndex_surfDispIncr>();
    
    constexpr MeshFieldIndex mfIndex_metric = MeshF_TanLatVertexRotatedOverRadius;
    auto meshField_metric = getMeshField<mfIndex_metric>();

    constexpr MeshFieldIndex mfIndex_vel = MeshF_Vel;
    auto meshFieldVel = getMeshField<mfIndex_vel>();

    Kokkos::parallel_for("update_velocity", numVtx, KOKKOS_LAMBDA(int vtx) {
      auto u =  meshFieldVel(vtx, 0);
      auto v =  meshFieldVel(vtx, 1);
      auto Del = elasticTimeStep * u * meshField_metric(vtx, 0);
      meshFieldVel(vtx, 0) =  cos(Del) * u + sin(Del) * v;
      meshFieldVel(vtx, 1) = -sin(Del) * u + cos(Del) * v;
    });

    Kokkos::parallel_for("calcVelIncr", numVtx, KOKKOS_LAMBDA(int vtx){
      meshField_surfDispIncr(vtx, 0) = meshField_surfDispIncr(vtx, 0) + meshFieldVel(vtx, 0);
      meshField_surfDispIncr(vtx, 1) = meshField_surfDispIncr(vtx, 1) + meshFieldVel(vtx, 1);
    });
  }

  void Mesh::applyFreeSlipBC(){
    Kokkos::Timer timer;
    auto nVerticesOwned= getNumVerticesOwned();
    auto meshFieldVel = getMeshField<MeshF_Vel>();
    auto uBdryVtxNormal = getMeshField<MeshF_uBdryVtxNormal>();
    auto vBdryVtxNormal = getMeshField<MeshF_vBdryVtxNormal>();
    auto interiorVertex = getMeshField<MeshF_InteriorVertex>();
    Kokkos::parallel_for("set free slip BC", nVerticesOwned, KOKKOS_LAMBDA(const int iVtx){
      if(interiorVertex(iVtx)) return;

      // edge 1 outward normal velocity
      auto normalVelocity1 = meshFieldVel(iVtx, 0) * uBdryVtxNormal(iVtx, 0) + meshFieldVel(iVtx, 1) * vBdryVtxNormal(iVtx, 0);

      // edge 2 outward normal velocity
      auto normalVelocity2 = meshFieldVel(iVtx, 0) * uBdryVtxNormal(iVtx, 1) + meshFieldVel(iVtx, 1) * vBdryVtxNormal(iVtx, 1);

      // remove edge-1 normal component
      auto uTry1 = meshFieldVel(iVtx, 0) - normalVelocity1 * uBdryVtxNormal(iVtx, 0);
      auto vTry1 = meshFieldVel(iVtx, 1) - normalVelocity1 * vBdryVtxNormal(iVtx, 0);

      // remove edge-2 normal component
      auto uTry2 = meshFieldVel(iVtx, 0) - normalVelocity2 * uBdryVtxNormal(iVtx, 1);
      auto vTry2 = meshFieldVel(iVtx, 1) - normalVelocity2 * vBdryVtxNormal(iVtx, 1);

      auto normalVelocity2WithTry1 = uTry1 * uBdryVtxNormal(iVtx, 1) + vTry1 * vBdryVtxNormal(iVtx, 1);

      auto normalVelocity1WithTry2 = uTry2 * uBdryVtxNormal(iVtx, 0) + vTry2 * vBdryVtxNormal(iVtx, 0);

      // already satisfies free-slip constraint
      if (normalVelocity1 <= 0.0 && normalVelocity2 <= 0.0){
        return;
      }
      else if (normalVelocity1 > 0.0 && normalVelocity2WithTry1 <= 0.0){
        meshFieldVel(iVtx,0) = uTry1;
        meshFieldVel(iVtx,1) = vTry1;
      }
      else if (normalVelocity2 > 0.0 && normalVelocity1WithTry2 <= 0.0){
        meshFieldVel(iVtx,0) = uTry2;
        meshFieldVel(iVtx,1) = vTry2;
      }
      else{
        meshFieldVel(iVtx,0) = 0.0;
        meshFieldVel(iVtx,1) = 0.0;
      }
    });
    pumipic::RecordTime("free slip BC", timer.seconds());
  }

} // namespace polyMPO
