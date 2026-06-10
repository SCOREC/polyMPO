#include <Kokkos_Core.hpp>
#include "pmpo_utils.hpp"
#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"
#include "pmpo_const_relation.hpp"

namespace polyMPO{

void printVTP_mesh(MPMesh& mpMesh, int printVTPIndex=-1);

void MPMesh::calculateStrain(){
  auto MPsPosition = p_MPs->getPositions();
  auto MPsBasis = p_MPs->getData<MPF_Basis_Vals>();
  auto MPsBasisGrads = p_MPs->getData<MPF_Basis_Grad_Vals>();
  auto MPsAppID = p_MPs->getData<MPF_MP_APP_ID>();
  auto MPsStrainRate = p_MPs->getData<MPF_Strain_Rate>();
  //Mesh Fields
  auto tanLatVertexRotatedOverRadius = p_mesh->getMeshField<MeshF_TanLatVertexRotatedOverRadius>();
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto velField = p_mesh->getMeshField<MeshF_Vel>();
  auto solveStress = p_mesh->getMeshField<polyMPO::MeshF_SolveStress>();

  auto setMPStrainRate = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){

      if(solveStress(elm)==0){
        MPsStrainRate(mp, 0) =  0.0;
        MPsStrainRate(mp, 1) =  0.0;
        MPsStrainRate(mp, 2) =  0.0;
        return;
      }

      int numVtx = elm2VtxConn(elm,0);

      double v11 = 0.0;
      double v12 = 0.0;
      double v21 = 0.0;
      double v22 = 0.0;
      double uTanOverR = 0.0;
      double vTanOverR = 0.0;
 
      for (int i = 0; i < numVtx; i++){
        int iVertex = elm2VtxConn(elm, i+1)-1;
        v11 = v11 + MPsBasisGrads(mp, i*2 + 0) * velField(iVertex, 0);
        v12 = v12 + MPsBasisGrads(mp, i*2 + 1) * velField(iVertex, 0);
        v21 = v21 + MPsBasisGrads(mp, i*2 + 0) * velField(iVertex, 1);
        v22 = v22 + MPsBasisGrads(mp, i*2 + 1) * velField(iVertex, 1);
        uTanOverR = uTanOverR + MPsBasis(mp, i) * tanLatVertexRotatedOverRadius(iVertex, 0) * velField(iVertex, 0);
        vTanOverR = vTanOverR + MPsBasis(mp, i) * tanLatVertexRotatedOverRadius(iVertex, 0) * velField(iVertex, 1);
      }

      MPsStrainRate(mp, 0) =  v11 - vTanOverR;
      MPsStrainRate(mp, 1) =  v22;
      MPsStrainRate(mp, 2) =  0.5*(v12 + v21 + uTanOverR);
    }
  };
  p_MPs->parallel_for(setMPStrainRate, "setMPStrainRate");
}

void MPMesh::calculateStress(const int constitutive_relation){
  //MeshFields  
  auto solveStress = p_mesh->getMeshField<polyMPO::MeshF_SolveStress>();
  auto elasticTimeStep = p_mesh->getElasticTimeStep();
  auto dynamicTimeStep = p_mesh->getDynamicTimeStep();
  auto dampingTimescale = polyMPO::dampingTimescaleParameter * dynamicTimeStep;
  //MPFields
  auto MPsAppID       = p_MPs->getData<MPF_MP_APP_ID>();
  auto MPsStrainRate  = p_MPs->getData<MPF_Strain_Rate>();
  auto MPsStress      = p_MPs->getData<MPF_Stress>();
  auto MPsArea        = p_MPs->getData<polyMPO::MPF_Area>();
  auto MPsIcePressure = p_MPs->getData<polyMPO::MPF_IcePressure>();
  auto MPsRepPressure = p_MPs->getData<polyMPO::MPF_ReplacementPressure>();

  auto setMPStress = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
    if(mask){
      Vec3d strain_rate (MPsStrainRate(mp, 0), MPsStrainRate(mp, 1), MPsStrainRate(mp, 2));
      Vec3d stress(MPsStress(mp, 0), MPsStress(mp, 1), MPsStress(mp, 2));

      if (constitutive_relation == 1)
        constitutive_evp(strain_rate, stress, MPsIcePressure(mp,0), MPsRepPressure(mp,0), MPsArea(mp,0), elasticTimeStep, dampingTimescale);
      else if(constitutive_relation == 3)
        constitutive_linear(strain_rate, stress);
      for (int m=0 ; m<3; m++)
        MPsStress(mp, m) = stress[m]*solveStress(elm);
    }
  };
  p_MPs->parallel_for(setMPStress, "setMPStress");
}

void MPMesh::calculateStressDivergence(){

  Kokkos::Timer timer;
  int self, numProcsTot;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  MPI_Comm_size(comm, &numProcsTot);

  //Mesh Information
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  int numVtxOwned = p_mesh->getNumVerticesOwned();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  int numVertices = p_mesh->getNumVertices();
  auto tanLatVertexRotatedOverRadius = p_mesh->getMeshField<MeshF_TanLatVertexRotatedOverRadius>();
  auto interiorVertex = p_mesh->getMeshField<MeshF_InteriorVertex>();

  //Material Points
  auto MPsAppID  = p_MPs->getData<MPF_MP_APP_ID>();
  auto weight = p_MPs->getData<MPF_Basis_Vals>();
  auto weight_grads = p_MPs->getData<MPF_Basis_Grad_Vals>(); 
  auto mpPos = p_MPs->getData<MPF_Cur_Pos_XYZ>();
  auto MPsStress = p_MPs->getData<MPF_Stress>();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>(); 

  auto VtxCoeffs_new   = this->precomputedVtxCoeffs_new;
  auto vtxMatrixMass_l = this->vtxMatrixMass;
  auto nearAnEdge_l    = this->nearAnEdge;

  //Earth Radius
  double radius = 1.0;
  if(p_mesh->getGeomType() == geom_spherical_surf)
    radius=p_mesh->getSphereRadius();

  auto stress_divUV = p_mesh->getMeshField<MeshF_StressDivergence>();
  Kokkos::deep_copy(stress_divUV, 0.0);

  //Assemble fields for Stress Divergence
  auto stress_div = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1;
        
        double ramp = nearAnEdge_l(vID);
        double invM = 1.0/vtxMatrixMass_l(vID);
        invM = vtxMatrixMass_l(vID) >1e-4 ? invM : 0; 

        double w_vtx=weight(mp,i); 
        double CoordDiffs[vec4d_nEntries] = {1, (-vtxCoords(vID,0) + mpPositions(mp,0))/radius,
                                                (-vtxCoords(vID,1) + mpPositions(mp,1))/radius,
                                                (-vtxCoords(vID,2) + mpPositions(mp,2))/radius};

        auto factor = ramp * w_vtx * (VtxCoeffs_new(vID,0, 0) + VtxCoeffs_new(vID,0, 1)*CoordDiffs[1]  +
                                                                VtxCoeffs_new(vID,0, 2)*CoordDiffs[2]  +
                                                                VtxCoeffs_new(vID,0, 3)*CoordDiffs[3]) +
                                                                (1.0 - ramp) * invM * w_vtx;

        factor = factor * tanLatVertexRotatedOverRadius(vID, 0);
      
        auto factor1 = ramp * (w_vtx/radius) * (VtxCoeffs_new(vID, 1, 0) + VtxCoeffs_new(vID, 1, 1)*CoordDiffs[1]  +
                                                                           VtxCoeffs_new(vID, 1, 2)*CoordDiffs[2]  +
                                                                           VtxCoeffs_new(vID, 1, 3)*CoordDiffs[3]) -
                                                                           (1.0 - ramp) * invM * weight_grads(mp, i*2+0);

        auto factor2 = ramp * (w_vtx/radius) * (VtxCoeffs_new(vID, 2, 0) + VtxCoeffs_new(vID, 2, 1)*CoordDiffs[1]  +
                                                                           VtxCoeffs_new(vID, 2, 2)*CoordDiffs[2]  +
                                                                           VtxCoeffs_new(vID, 2, 3)*CoordDiffs[3]) -
                                                                           (1.0 - ramp) * invM * weight_grads(mp, i*2+1);

        Kokkos::atomic_add(&stress_divUV(vID, 0), factor1 * MPsStress(mp, 0) + factor2 * MPsStress(mp, 2) -
                                                  2 * factor * MPsStress(mp, 2));

        Kokkos::atomic_add(&stress_divUV(vID, 1), factor2 * MPsStress(mp, 1) + factor1 * MPsStress(mp, 2) +
                                                  factor * (MPsStress(mp, 0)-MPsStress(mp, 1)));
      }
    }
  };
  p_MPs->parallel_for(stress_div, " stress_div_assembly");
  Kokkos::fence();
  pumipic::RecordTime("Stress_Divergence_Reconstruction" + std::to_string(self), timer.seconds()); 

  timer.reset();
  if(numProcsTot>1){ 
    //Takes contribution of halo vertices and adds it in owner procs
    communicate_and_take_halo_contributions1(stress_divUV, numVertices, 2, 0, 0);
    //Transfer the correct values at owned vertices to halo vertices
    //communicate_and_take_halo_contributions(stress_divUV, numVertices, 2, 1, 1);
  }
  Kokkos::fence();
  pumipic::RecordTime("Stress_Divergence Communication" + std::to_string(self), timer.seconds());  
}

void MPMesh::calcBasis() {
  assert(p_mesh->getGeomType() == geom_spherical_surf);

  auto MPsPosition = p_MPs->getPositions();
  auto MPsBasis = p_MPs->getData<MPF_Basis_Vals>();
  auto MPsBasisGrads = p_MPs->getData<MPF_Basis_Grad_Vals>();
  auto MPsAppID = p_MPs->getData<MPF_MP_APP_ID>();

  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto vtxCoords = p_mesh->getMeshField<MeshF_VtxCoords>();
  double radius = 1.0;
  if(p_mesh->getGeomType() == geom_spherical_surf)
    radius=p_mesh->getSphereRadius();

  //For Gnomonic Projection
  auto gnomProjVtx = p_mesh->getMeshField<polyMPO::MeshF_VtxGnomProj>();
  auto gnomProjElmCenter = p_mesh->getMeshField<polyMPO::MeshF_ElmCenterGnomProj>();

  bool isRotated = p_mesh->getRotatedFlag();

  auto calcbasis = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int numVtx = elm2VtxConn(elm,0);
      Vec3d position3d(MPsPosition(mp, 0),MPsPosition(mp, 1),MPsPosition(mp, 2));
      if(isRotated){
        position3d[0] = -MPsPosition(mp, 2);
        position3d[2] = MPsPosition(mp, 0);
      }

      double mpProjX, mpProjY;
      auto gnomProjElmCenter_sub = Kokkos::subview(gnomProjElmCenter, elm, Kokkos::ALL);
      computeGnomonicProjectionAtPoint(position3d, gnomProjElmCenter_sub, mpProjX, mpProjY);
      auto gnom_vtx_subview = Kokkos::subview(gnomProjVtx, elm, Kokkos::ALL, Kokkos::ALL); 

      double basisByArea[maxVtxsPerElm] = {0.0};
      initArray(basisByArea,maxVtxsPerElm, 0.0);
      double gradBasisByArea[2*maxVtxsPerElm] = {0.0};
      initArray(gradBasisByArea,maxVtxsPerElm, 0.0);

      wachpress_weights_grads_2D(numVtx, gnom_vtx_subview, mpProjX, mpProjY, radius, basisByArea, gradBasisByArea);

      for(int i=0; i<= numVtx; i++){
        MPsBasis(mp, i) = basisByArea[i];
        MPsBasisGrads(mp, i*2+0) = gradBasisByArea[i*2 + 0];
        MPsBasisGrads(mp, i*2+1) = gradBasisByArea[i*2 + 1];
      }

      //Old method where basis functions calculated using 3D Area
      /*
      Vec3d v3d[maxVtxsPerElm+1];
      int numVtx = elm2VtxConn(elm,0);
      for(int i = 1; i<=numVtx; i++){
        v3d[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
        v3d[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
        v3d[i-1][2] = vtxCoords(elm2VtxConn(elm,i)-1,2);
      }
      v3d[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
      v3d[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
      v3d[numVtx][2] = vtxCoords(elm2VtxConn(elm,1)-1,2); 
      getBasisByAreaGblFormSpherical(position3d, numVtx, v3d, radius, basisByArea);
      */
    }
  };
  p_MPs->parallel_for(calcbasis, "calcbasis");
}

void MPMesh::CVTTrackingElmCenterBased(const int printVTPIndex){
  Kokkos::Timer timer;
  int numVtxs = p_mesh->getNumVertices();
  int numElms = p_mesh->getNumElements();
  auto numMPs = p_MPs->getCount();

  const auto elmCenter = p_mesh->getMeshField<polyMPO::MeshF_ElmCenterXYZ>();

  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto elm2ElmConn = p_mesh->getElm2ElmConn();

  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
  auto mpTgtPos = p_MPs->getData<MPF_Tgt_Pos_XYZ>();
  auto MPs2Elm = p_MPs->getData<MPF_Tgt_Elm_ID>();
  auto MPs2Proc = p_MPs->getData<MPF_Tgt_Proc_ID>();
  auto elm2Process = p_mesh->getElm2Process();
  auto elm2global = p_mesh->getElmGlobal();

  if(printVTPIndex>=0) {
    printVTP_mesh(printVTPIndex);
  }

  Vec3dView history("positionHistory",numMPs);
  Vec3dView resultLeft("positionResult",numMPs);
  Vec3dView resultRight("positionResult",numMPs);
  Vec3dView mpTgtPosArray("positionTarget",numMPs);

  auto CVTElmCalc = PS_LAMBDA(const int& elm, const int& mp, const int&mask){
    Vec3d MP(mpPositions(mp,0),mpPositions(mp,1),mpPositions(mp,2));
    if(mask){
      int iElm = elm;
      Vec3d MPnew(mpTgtPos(mp,0),mpTgtPos(mp,1),mpTgtPos(mp,2));
      Vec3d dx = MPnew-MP;
      while(true){
        int numConnElms = elm2ElmConn(iElm,0);
                
        Vec3d center(elmCenter(iElm, 0), elmCenter(iElm, 1), elmCenter(iElm, 2));
        Vec3d delta = MPnew - center;

        double minDistSq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
        int closestElm = -1;
        //go through all the connected elm, calc distance
        for(int i=1; i<=numConnElms; i++){
          int elmID = elm2ElmConn(iElm,i)-1;
          if (elmID >= numElms)
            continue; 
          //New delta
          Vec3d center(elmCenter(elmID, 0), elmCenter(elmID, 1), elmCenter(elmID, 2));
          delta = MPnew - center;

          double neighborDistSq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
          if(neighborDistSq < minDistSq){
            closestElm = elmID;
            minDistSq = neighborDistSq;
          }
        }

        if(closestElm<0){
          MPs2Elm(mp) = iElm;
          if (elm2Process.size() > 0)
            MPs2Proc(mp) = elm2Process(iElm);
          break;
        }
        else{
          iElm = closestElm;
        }
      }
        
      if(printVTPIndex>=0 && numMPs>0){
        double d1 = dx[0];
        double d2 = dx[2];
        double d3 = dx[3];
        double m1 = MP[0];
        double m2 = MP[1];
        double m3 = MP[2];
        Vec3d MParrow = MP + dx*0.7;
        Vec3d r = MPnew * (1.0/MPnew.magnitude());
        Vec3d shift = dx.cross(r) * ((1.0-0.7)*dx.magnitude()/(dx.cross(r)).magnitude());
        Vec3d MPLeft = MParrow + shift;
        Vec3d MPRight = MParrow - shift;
        history(mp) = MP;
        resultLeft(mp) = MPLeft;
        resultRight(mp) = MPRight;
        mpTgtPosArray(mp) = MPnew;
      }
    }
  };
  p_MPs->parallel_for(CVTElmCalc,"CVTTrackingElmCenterBasedCalc");

  if(printVTPIndex>=0){
    writeMPTrackingVTP(printVTPIndex, numMPs, history, resultLeft, resultRight, mpTgtPosArray);
  }
  pumipic::RecordTime("PolyMPO_CVTTrackingElmCenterBased", timer.seconds());
}

void MPMesh::reconstructSlices() {
  if (reconstructSlice.size() == 0) return;
  Kokkos::Timer timer;
  for (auto const& [index, reconstruct] : reconstructSlice) {
    if (reconstruct) reconstruct();
  }
  reconstructSlice.clear();
  pumipic::RecordTime("PolyMPO_Reconstruct", timer.seconds());
}

bool getAnyIsMigrating(MaterialPoints* p_MPs, bool isMigrating) {
  Kokkos::Timer timer;
  MPI_Comm comm = p_MPs->getMPIComm();
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);
  int comm_size;
  MPI_Comm_size(comm, &comm_size);

  bool anyIsMigrating = false;
  MPI_Allreduce(&isMigrating, &anyIsMigrating, 1, MPI_C_BOOL, MPI_LOR, comm);
  pumipic::RecordTime("PolyMPO_getAnyIsMigrating", timer.seconds());
  return anyIsMigrating;
}

void MPMesh::push_ahead(){
  Kokkos::Timer timer;
  //Latitude Longitude increment at mesh vertices
  p_mesh->computeRotLatLonIncr();   

  //Interpolates latitude longitude, mesh velocity increments to MPs
  //calcBasis();
  sphericalInterpolation<MeshF_RotLatLonIncr>(*this);
  sphericalInterpolation<MeshF_OnSurfVeloIncr>(*this);
  sphericalInterpolation2Fields(*this);
  //Push the MPs
  p_MPs->updateRotLatLonAndXYZ2Tgt(p_mesh->getSphereRadius(), p_mesh->getRotatedFlag());
  pumipic::RecordTime("PolyMPO_interpolateAndPush", timer.seconds());
}

bool MPMesh::push1P(){
  Kokkos::Timer timer;
  //Given target location find the new element or the last element in a partioned mesh
  //and the process it belongs to so that migration can be checked
  CVTTrackingElmCenterBased(); 
  //From the above two inputs find if any particle needs to be migrated
  bool anyIsMigrating = getAnyIsMigrating(p_MPs, p_MPs->check_migrate());
  pumipic::RecordTime("PolyMPO_trackAndCheckMigrate", timer.seconds());
  return anyIsMigrating;
}

void MPMesh::push_swap(){
  //current becomes target, target becomes -1
  p_MPs->updateMPElmID();
}

void MPMesh::push_swap_pos(){
  //current becomes target, target becomes -1
  //Making read for next push_ahead  
  p_MPs->updateMPSlice<MPF_Cur_Pos_XYZ, MPF_Tgt_Pos_XYZ>();
  p_MPs->updateMPSlice<MPF_Cur_Pos_Rot_Lat_Lon, MPF_Tgt_Pos_Rot_Lat_Lon>();
}

void MPMesh::push(){  
  Kokkos::Timer timer;

  p_mesh->computeRotLatLonIncr();

  sphericalInterpolation<MeshF_RotLatLonIncr>(*this);

  p_MPs->updateRotLatLonAndXYZ2Tgt(p_mesh->getSphereRadius(), p_mesh->getRotatedFlag());

  auto elm2Process = p_mesh->getElm2Process();

  bool anyIsMigrating = false;
  do {
    CVTTrackingElmCenterBased(); // move to Tgt_XYZ
    p_MPs->updateMPSlice<MPF_Cur_Pos_XYZ, MPF_Tgt_Pos_XYZ>(); // Tgt_XYZ becomes Cur_XYZ
    p_MPs->updateMPSlice<MPF_Cur_Pos_Rot_Lat_Lon, MPF_Tgt_Pos_Rot_Lat_Lon>(); // Tgt becomes Cur
    
    bool anyIsMigrating = getAnyIsMigrating(p_MPs, p_MPs->check_migrate());

    if(anyIsMigrating)
      p_MPs->migrate();
    else
      p_MPs->rebuild();

    p_MPs->updateMPElmID(); //update mpElm IDs slices
    reconstructSlices(); 
  }
  while (anyIsMigrating);

  pumipic::RecordTime("PolyMPO_push", timer.seconds());
}

//Writing vtp Mesh and MP history
void MPMesh::printVTP_mesh(int printVTPIndex){
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  auto elm2VtxConn = p_mesh->getElm2VtxConn();

  auto MPsPosition = p_MPs->getPositions();

  char* fileOutput = (char *)malloc(sizeof(char) * 256); 
  sprintf(fileOutput,"polyMPO_MPMesh_mesh_%d.vtp", printVTPIndex);
  FILE * pFile = fopen(fileOutput,"w");
  free(fileOutput);

  auto h_vtxCoords = Kokkos::create_mirror_view(vtxCoords);
  IntVtx2ElmView::HostMirror h_elm2VtxConn = Kokkos::create_mirror_view(elm2VtxConn);
  const int nCells = p_mesh->getNumElements();
  const int nVertices = p_mesh->getNumVertices();
  Kokkos::deep_copy(h_vtxCoords,vtxCoords);
  Kokkos::deep_copy(h_elm2VtxConn,elm2VtxConn);
  fprintf(pFile, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%d\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",nVertices,nCells);
  for(int i=0; i<nVertices; i++){
    fprintf(pFile, "          %f %f %f\n",h_vtxCoords(i,0),h_vtxCoords(i,1),h_vtxCoords(i,2));
  }
  fprintf(pFile, "        </DataArray>\n      </Points>\n      <Polys>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0; i<nCells; i++){
    fprintf(pFile, "          ");
    for(int j=0; j< h_elm2VtxConn(i,0); j++){
      fprintf(pFile, "%d ", h_elm2VtxConn(i,j+1)-1);
    } 
    fprintf(pFile, "\n");
  }
  fprintf(pFile, "        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    
  int count = 0;
  for(int i=0;i<nCells; i++){
    count += h_elm2VtxConn(i,0);
    fprintf(pFile, "          %d\n",count);
  }
  fprintf(pFile, "        </DataArray>\n      </Polys>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
  fclose(pFile);
}

void MPMesh::writeMPTrackingVTP(int printVTPIndex, int numMPs, const Vec3dView& history, const Vec3dView& resultLeft,
                                const Vec3dView& resultRight, const Vec3dView& mpTgtPosArray){

  Vec3dView::HostMirror h_history = Kokkos::create_mirror_view(history);
  Vec3dView::HostMirror h_resultLeft = Kokkos::create_mirror_view(resultLeft);
  Vec3dView::HostMirror h_resultRight = Kokkos::create_mirror_view(resultRight);
  Vec3dView::HostMirror h_mpTgtPos = Kokkos::create_mirror_view(mpTgtPosArray);

  Kokkos::deep_copy(h_history, history);
  Kokkos::deep_copy(h_resultLeft, resultLeft);
  Kokkos::deep_copy(h_resultRight, resultRight);
  Kokkos::deep_copy(h_mpTgtPos, mpTgtPosArray);

  char* fileOutput = (char *)malloc(sizeof(char) * 256); 
  sprintf(fileOutput, "polyMPOCVTTrackingElmCenter_MPtracks_%d.vtp", printVTPIndex);
  FILE * pFile = fopen(fileOutput,"w");
  free(fileOutput);   

  fprintf(pFile, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"%d\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",numMPs*4,numMPs*2); 
  for(int i=0; i<numMPs; i++){
    fprintf(pFile,"          %f %f %f\n          %f %f %f\n          %f %f %f\n          %f %f %f\n",
      h_history(i)[0],h_history(i)[1],h_history(i)[2],
      h_mpTgtPos(i)[0],h_mpTgtPos(i)[1],h_mpTgtPos(i)[2],
      h_resultLeft(i)[0],h_resultLeft(i)[1],h_resultLeft(i)[2],
      h_resultRight(i)[0],h_resultRight(i)[1],h_resultRight(i)[2]);
  }
  fprintf(pFile,"        </DataArray>\n      </Points>\n      <Lines>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"); 
  for(int i=0; i<numMPs*4; i+=4){
    fprintf(pFile,"          %d %d\n          %d %d %d\n",i,i+1,i+2,i+1,i+3);
  }
  fprintf(pFile,"        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  for(int i=0; i<numMPs*5; i+=5){
    fprintf(pFile,"          %d\n          %d\n",i+2,i+5);
  }
  fprintf(pFile,"        </DataArray>\n      </Lines>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
  fclose(pFile);

}

//These are not used currently
void MPMesh::CVTTrackingEdgeCenterBased(Vec2dView dx){
    int numElms = p_mesh->getNumElements();

    auto elm2VtxConn = p_mesh->getElm2VtxConn();
    auto elm2ElmConn = p_mesh->getElm2ElmConn();
    auto MPs2Elm = p_MPs->getData<MPF_Tgt_Elm_ID>();
    const auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
    Kokkos::View<Vec2d*[maxVtxsPerElm]> edgeCenters("EdgeCenters",numElms);

    Kokkos::parallel_for("calcEdgeCenter", numElms, KOKKOS_LAMBDA(const int elm){  
        int numVtx = elm2VtxConn(elm,0);
        int v[maxVtxsPerElm];
        for(int i=0; i< numVtx; i++)
            v[i] = elm2VtxConn(elm,i+1)-1;
        for(int i=0; i< numVtx; i++){
            int idx_ip1 = (i+1)%numVtx;
            Vec2d v_i(vtxCoords(v[i],0),vtxCoords(v[i],1));
            Vec2d v_ip1(vtxCoords(v[idx_ip1],0),vtxCoords(v[idx_ip1],1));
            edgeCenters(elm,i) = (v_ip1 + v_i)*0.5;
        }
    });

    auto CVTEdgeTracking = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
        Vec2d MP(mpPositions(mp,0),mpPositions(mp,1));//XXX:the input is XYZ, but we only support 2d vector
        if(mask){
            Vec2d MPnew = MP + dx(mp);
            int iElm = elm;
            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                //calc dist square from each edge center to MPnew
                //calc dot products to check inside or not
                int edgeIndex = -1;
                double minDistSq = DBL_MAX;
                for(int i=0; i< numVtx; i++){
                    Vec2d edgeCenter = edgeCenters(iElm,i);
                    Vec2d delta = MPnew - edgeCenter;
                    double currentDistSq = delta[0]*delta[0] + delta[1]*delta[1];
                    double dotProduct = dx(mp).dot(delta);
                    if(dotProduct <=0){
                        edgeIndex = -1;
                        break;
                    }
                    if(currentDistSq < minDistSq){
                        edgeIndex = i+1;
                        minDistSq = currentDistSq;
                    }
                }
                if(edgeIndex <0){
                    //we get to the final elm
                    MPs2Elm(mp) = iElm; 
                    mpPositions(mp,0) = MPnew[0];
                    mpPositions(mp,1) = MPnew[1];
                    mpPositions(mp,2) = 0.0; //XXX:we only have 2d vector
                    break;
                }else{
                    //update the iELm and do the loop again
                    iElm = elm2ElmConn(iElm,edgeIndex);
                }
            } 
        }
    };
    p_MPs->parallel_for(CVTEdgeTracking,"CVTTrackingEdgeCenterBasedCalc");
}

void MPMesh::T2LTracking(Vec2dView dx){
    const auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>(); 
    auto elm2VtxConn = p_mesh->getElm2VtxConn();
    auto elm2ElmConn = p_mesh->getElm2ElmConn();

    auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
    auto MPs2Elm = p_MPs->getData<MPF_Tgt_Elm_ID>();
    auto mpStatus = p_MPs->getData<MPF_Status>();

    auto T2LCalc = PS_LAMBDA(const int& elm, const int& mp, const int&mask){
        Vec2d MP(mpPositions(mp,0),mpPositions(mp,1));//XXX:the input is XYZ, but we only support 2d vector
        if(mask){
            int iElm = elm;
            Vec2d MPnew = MP + dx(mp);    
            
            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                bool goToNeighbour = false;
                //seperate the elm2Vtx
                int v[maxVtxsPerElm];
                for(int i=0; i< numVtx; i++)
                    v[i] = elm2VtxConn(iElm,i+1)-1;
                //get edges and perpendiculardx
                Vec2d e[maxVtxsPerElm];
                double pdx[maxVtxsPerElm];                    
                for(int i=0; i< numVtx; i++){
                    int idx_ip1 = (i+1)%numVtx;
                    Vec2d v_i(vtxCoords(v[i],0),vtxCoords(v[i],1));
                    Vec2d v_ip1(vtxCoords(v[idx_ip1],0),vtxCoords(v[idx_ip1],1));
                    e[i] = v_ip1 - v_i;
                    pdx[i] = (v_i - MP).cross(dx(mp));
                }
                
                for(int i=0; i<numVtx; i++){
                    int ip1 = (i+1)%numVtx;
                    //pdx*pdx<0 and edge is acrossed 
                    if(pdx[i]*pdx[ip1] <0 && e[i].cross(Vec2d(MPnew[0]-vtxCoords(v[i],0),
                                                              MPnew[1]-vtxCoords(v[i],1)))<0){
                        //go to the next elm
                        iElm = elm2ElmConn(iElm,i+1);
                        goToNeighbour = true;
                        if(iElm <0){
                            mpStatus(mp) = 0;                  
                            MPs2Elm(mp) = -1;
                            goToNeighbour = false;
                        }
                    }
                }
                //if goes to the other 
                if(goToNeighbour)
                    continue; 
                //otherwise we do the update and end the loop
                MPs2Elm(mp) = iElm;
                mpPositions(mp,0) = MPnew[0];
                mpPositions(mp,1) = MPnew[1];
                mpPositions(mp,2) = 0.0; //XXX:we only have 2d vector
                break;
            }
        }
    }; 
    p_MPs->parallel_for(T2LCalc,"T2lTrackingCalc");
}

}
