#ifndef POLYMPO_ASSEMBLY_H
#define POLYMPO_ASSEMBLY_H

#include "pmpo_wachspressBasis.hpp"

namespace polyMPO{

DoubleView MPMesh::assemblyV0(){
  int numVtxs = p_mesh->getNumVertices();
  auto elm2VtxConn = p_mesh->getElm2VtxConn();

  DoubleView vField("vField2",numVtxs);
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>(); //get the array of MP coordinates/positions
  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
        double distance = mpPositions(mp,0) + mpPositions(mp,1) + mpPositions(mp,2);
        Kokkos::atomic_add(&vField(vID),distance);
      }
    }
  };
  p_MPs->parallel_for(assemble, "assembly");
  return vField;
}

template <MeshFieldIndex meshFieldIndex>
void MPMesh::assemblyVtx0(){
  Kokkos::Timer timer;

  constexpr MaterialPointSlice mpfIndex = meshFieldIndexToMPSlice<meshFieldIndex>;
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto mpData = p_MPs->getData<mpfIndex>();
  const int numEntries = mpSliceToNumEntries<mpfIndex>();

  int numVtx = p_mesh->getNumVertices();
  p_mesh->fillMeshField<meshFieldIndex>(numVtx, numEntries, 0.0);
  auto meshField = p_mesh->getMeshField<meshFieldIndex>();
  auto weight = p_MPs->getData<MPF_Basis_Vals>();

  const double zero = 0.0;
  Kokkos::View<double*> sumWeights("sumWeights", numVtx);
  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
        double fieldComponentVal;
        Kokkos::atomic_add(&sumWeights(vID), weight(mp, 0));
        for(int j=0;j<numEntries;j++){
          fieldComponentVal = mpData(mp,j) * weight(mp, 0);
          Kokkos::atomic_add(&meshField(vID,j),fieldComponentVal);
        }
      }
    }
  };
  p_MPs->parallel_for(assemble, "assembly");
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},{numVtx, numEntries});
  Kokkos::parallel_for("assembly average", policy, KOKKOS_LAMBDA(const int vtx, const int entry){
    if (sumWeights(vtx) != zero) 
      meshField(vtx, entry) /= sumWeights(vtx);
  });
  pumipic::RecordTime("PolyMPO_Reconstruct_Vtx0", timer.seconds());
}

template <MeshFieldIndex meshFieldIndex>
void MPMesh::assemblyElm0() {
  Kokkos::Timer timer;
  constexpr MaterialPointSlice mpfIndex = meshFieldIndexToMPSlice<meshFieldIndex>;
  auto mpData = p_MPs->getData<mpfIndex>();
  const int numEntries = mpSliceToNumEntries<mpfIndex>();

  int numElms = p_mesh->getNumElements();
  p_mesh->fillMeshField<meshFieldIndex>(numElms, numEntries, 0.0);
  auto meshField = p_mesh->getMeshField<meshFieldIndex>();

  Kokkos::View<int*> mpsPerElm("mpsPerElm", numElms);
  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      Kokkos::atomic_add(&mpsPerElm(elm),1);
      for(int j=0;j<numEntries;j++){
        Kokkos::atomic_add(&meshField(elm,j), mpData(mp,j));
      }
    }
  };
  p_MPs->parallel_for(assemble, "assembly");
  
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},{numElms, numEntries});
  Kokkos::parallel_for("assembly average", policy, KOKKOS_LAMBDA(const int elm, const int entry){
    if (mpsPerElm(elm) > 0){
      meshField(elm, entry) /= mpsPerElm(elm);
    }
  });
  pumipic::RecordTime("PolyMPO_Reconstruct_Elm0", timer.seconds());
}

void MPMesh::reconstruct_coeff_full(){
  Kokkos::Timer timer;
  int self, numProcsTot;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  MPI_Comm_size(comm, &numProcsTot);
  
  static int coeff_count=0;
  if(!self) std::cout<<"===="<<__FUNCTION__<<" "<<coeff_count<<"===="<<std::endl;
  coeff_count++;
  //Mesh Information
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  int numVtx = p_mesh->getNumVertices();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  int numVertices = p_mesh->getNumVertices();
  //Dual Element Area for Regularization
  auto dual_triangle_area=p_mesh->getMeshField<MeshF_DualTriangleArea>();

  //Material Points
  calcBasis();

  auto weight = p_MPs->getData<MPF_Basis_Vals>();
  auto mpPos = p_MPs->getData<MPF_Cur_Pos_XYZ>();

  //Matrix for each vertex
  constexpr int numEntriesMatrix=10;
  Kokkos::View<double*[numEntriesMatrix]> vtxMatrices("VtxMatrices", p_mesh->getNumVertices());
  Kokkos::deep_copy(vtxMatrices, 0);

  //Earth Radius
  double radius = 1.0;
  if(p_mesh->getGeomType() == geom_spherical_surf)
    radius=p_mesh->getSphereRadius();

  //Assemble matrix for each vertex
  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
        double w_vtx=weight(mp,i);

        Kokkos::atomic_add(&vtxMatrices(vID,0), w_vtx);
        Kokkos::atomic_add(&vtxMatrices(vID,1), w_vtx*(-vtxCoords(vID,0)+mpPos(mp,0))/radius);
        Kokkos::atomic_add(&vtxMatrices(vID,2), w_vtx*(-vtxCoords(vID,1)+mpPos(mp,1))/radius);
        Kokkos::atomic_add(&vtxMatrices(vID,3), w_vtx*(-vtxCoords(vID,2)+mpPos(mp,2))/radius);
        Kokkos::atomic_add(&vtxMatrices(vID,4), w_vtx*(-vtxCoords(vID,0)+mpPos(mp,0))*(-vtxCoords(vID,0)+mpPos(mp,0))/(radius*radius));
        Kokkos::atomic_add(&vtxMatrices(vID,5), w_vtx*(-vtxCoords(vID,0)+mpPos(mp,0))*(-vtxCoords(vID,1)+mpPos(mp,1))/(radius*radius));
        Kokkos::atomic_add(&vtxMatrices(vID,6), w_vtx*(-vtxCoords(vID,0)+mpPos(mp,0))*(-vtxCoords(vID,2)+mpPos(mp,2))/(radius*radius));
        Kokkos::atomic_add(&vtxMatrices(vID,7), w_vtx*(-vtxCoords(vID,1)+mpPos(mp,1))*(-vtxCoords(vID,1)+mpPos(mp,1))/(radius*radius));
        Kokkos::atomic_add(&vtxMatrices(vID,8), w_vtx*(-vtxCoords(vID,1)+mpPos(mp,1))*(-vtxCoords(vID,2)+mpPos(mp,2))/(radius*radius));
        Kokkos::atomic_add(&vtxMatrices(vID,9), w_vtx*(-vtxCoords(vID,2)+mpPos(mp,2))*(-vtxCoords(vID,2)+mpPos(mp,2))/(radius*radius));
      }
    }
  };
  p_MPs->parallel_for(assemble, "assembly");
  Kokkos::fence();
  pumipic::RecordTime("Assemble Matrix Per Process" + std::to_string(self), timer.seconds());
  //Mode 0 is Gather:  Halos Send to Owners
  //Mode 1 is Scatter: Owners Send to Halos
  //Op 0 is addition
  //Op 1 is replacement
  timer.reset();
  int mode = 0;
  int op = 0;
  if (numProcsTot >1){
    communicate_and_take_halo_contributions(vtxMatrices, numVertices, numEntriesMatrix, mode, op);
    mode=1; 
    op=1;
    communicate_and_take_halo_contributions(vtxMatrices, numVertices, numEntriesMatrix, mode, op);
  }
  pumipic::RecordTime("Communicate Matrix Values" + std::to_string(self), timer.seconds());
 
  //Stroe the 1st matrix element
  Kokkos::View<double*>vtxMatrixMass_l("vtxMass", numVertices);
  Kokkos::parallel_for("storeMatrixMass", numVertices, KOKKOS_LAMBDA(const int vtx){
    vtxMatrixMass_l(vtx) = vtxMatrices(vtx, 0);
  }); 
  this->vtxMatrixMass = vtxMatrixMass_l;

  invertMatrix(vtxMatrices, radius);
}

void MPMesh::invertMatrix(const Kokkos::View<double**>& vtxMatrices, const double& radius){
  
  int nVertices = p_mesh->getNumVertices();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  auto dual_triangle_area = p_mesh->getMeshField<MeshF_DualTriangleArea>();
  auto interiorVertex = p_mesh->getMeshField<MeshF_InteriorVertex>();
  bool isRotated = p_mesh->getRotatedFlag();

  double eps = 1e-7;
  double truncateFactor = 0.05;

  Kokkos::View<double*[3][vec4d_nEntries]> VtxCoeffs("VtxCoeffs", nVertices);
  Kokkos::deep_copy(VtxCoeffs, 0.0);
  Kokkos::View<double*> nearAnEdge_l("nearAnEdge_l", nVertices);
  
  Kokkos::parallel_for("invertMatrix", nVertices, KOKKOS_LAMBDA(const int vtx){
    if(vtxMatrices(vtx, 0) < eps)
      return;

    auto small = eps * vtxMatrices(vtx, 0) * dual_triangle_area(vtx, 0)/(radius*radius);
    auto truncate = truncateFactor * vtxMatrices(vtx, 0) * dual_triangle_area(vtx, 0)/(radius*radius);

    double X = vtxCoords(vtx, 0)/radius;
    double Y = vtxCoords(vtx, 1)/radius;
    double Z = vtxCoords(vtx, 2)/radius;
    if(isRotated){
      X = -vtxCoords(vtx, 2)/radius;
      Z = vtxCoords(vtx, 0)/radius;
    }

    auto cosLat = sqrt(pow(X, 2) +  pow(Y, 2));
    auto invCosLat = 1.0/cosLat;     
    auto vtx_area_sqrt = sqrt(dual_triangle_area(vtx,0)/(radius*radius));

    Vec3d v0 = { -Y * invCosLat,  -Z * X * invCosLat,  X / vtx_area_sqrt };
    Vec3d v1 = { X * invCosLat,  -Y * Z * invCosLat,  Y / vtx_area_sqrt };
    Vec3d v2 = { 0.0, cosLat, Z /vtx_area_sqrt };
    if(isRotated){
      v0 = {0.0, cosLat, Z / vtx_area_sqrt};
      v1 = {X * invCosLat, -Y * Z * invCosLat, Y / vtx_area_sqrt};
      v2 = {Y * invCosLat, X * Z *invCosLat, -X / vtx_area_sqrt};
    }
    Matrix3d rotateScaleM = {v0, v1, v2};

    double invM11 = 1.0 / vtxMatrices(vtx, 0);
    Matrix3d subM;
    subM(0, 0) = vtxMatrices(vtx, 4) - invM11 * vtxMatrices(vtx, 1) * vtxMatrices(vtx, 1);
    subM(0, 1) = vtxMatrices(vtx, 5) - invM11 * vtxMatrices(vtx, 1) * vtxMatrices(vtx, 2);
    subM(0, 2) = vtxMatrices(vtx, 6) - invM11 * vtxMatrices(vtx, 1) * vtxMatrices(vtx, 3);
    subM(1, 1) = vtxMatrices(vtx, 7) - invM11 * vtxMatrices(vtx, 2) * vtxMatrices(vtx, 2);
    subM(1, 2) = vtxMatrices(vtx, 8) - invM11 * vtxMatrices(vtx, 2) * vtxMatrices(vtx, 3);
    subM(2, 2) = vtxMatrices(vtx, 9) - invM11 * vtxMatrices(vtx, 3) * vtxMatrices(vtx, 3);
    subM(1, 0) = subM(0, 1);
    subM(2, 0) = subM(0, 2);
    subM(2, 1) = subM(1, 2);

    auto subM1 = (rotateScaleM.transpose())*(subM*rotateScaleM);

    Vec3d mVec = {vtxMatrices(vtx, 1), vtxMatrices(vtx, 2), vtxMatrices(vtx, 3)};
    auto blockC = rotateScaleM * mVec;

    auto trG = subM1(0, 0) + subM1(1, 1);
    if((trG < small) || (vtxMatrices(vtx, 0) < truncateFactor)){
      VtxCoeffs(vtx, 0, 0) = invM11;
      return;
    }

    auto diffTr = subM1(0, 0)-subM1(1, 1);
    auto delG = sqrt(4.0 * pow(subM1(0, 1), 2) + pow(diffTr, 2));
    if(delG<small){
      subM1(0, 0) = Kokkos::max(subM1(0, 0), truncate);
      subM1(1, 1) = Kokkos::max(subM1(1, 1), truncate);
      subM1(0, 1) = 0.0;
    }
    else{
      auto minEig = Kokkos::max(0.5 * (trG - delG), truncate);
      auto maxEig = Kokkos::max(0.5 * (trG + delG), truncate);
      auto trG = minEig + maxEig;
      auto diffEig = maxEig - minEig;
      diffTr = diffTr / delG;
      subM1(0, 0) = 0.5 * (trG + diffTr * diffEig);
      subM1(1, 1) = 0.5 * (trG - diffTr * diffEig);
      subM1(0, 1) = (1.0 / delG) * subM1(0, 1) * diffEig;
    }

    double denom = subM1(0, 0) * subM1(1, 1) - pow(subM1(0, 1), 2);
    double minZ2 = (subM1(1, 1) * pow(subM1(0, 2), 2) + subM1(0, 0) * pow(subM1(1, 2), 2) - 2.0 * subM1(0, 1) * subM1(0, 2) * subM1(1, 2))/denom;
    subM1(2, 2) = Kokkos::max(subM1(2, 2), abs(minZ2) +  truncate);

    double invM2D[6]={0.0};
    invM2D[0] = subM1(2, 2) * subM1(1, 1) - subM1(1, 2) * subM1(1, 2);
    invM2D[1] = subM1(0, 2) * subM1(1, 2) - subM1(2, 2) * subM1(0, 1);
    invM2D[2] = subM1(0, 1) * subM1(1, 2) - subM1(0, 2) * subM1(1, 1);
    invM2D[3] = subM1(2, 2) * subM1(0, 0) - subM1(0, 2) * subM1(0, 2);
    invM2D[4] = subM1(0, 1) * subM1(0, 2) - subM1(0, 0) * subM1(1, 2);
    invM2D[5] = subM1(0, 0) * subM1(1, 1) - subM1(0, 1) * subM1(0, 1);
    double det = subM1(0, 0) * invM2D[0] + subM1(0, 1) * invM2D[1] + subM1(0, 2) * invM2D[2];
    for (int i=0; i<6; i++) invM2D[i] = invM2D[i]/det;

    Vec3d iBlockC(0, 0, 0);
    iBlockC[0] = blockC[0] * invM2D[0] +  blockC[1] * invM2D[1] +  blockC[2] * invM2D[2];
    iBlockC[1] = blockC[0] * invM2D[1] +  blockC[1] * invM2D[3] +  blockC[2] * invM2D[4];
    iBlockC[2] = blockC[0] * invM2D[2] +  blockC[1] * invM2D[4] +  blockC[2] * invM2D[5];
    iBlockC = iBlockC*invM11; 

    VtxCoeffs(vtx, 0, 0) = invM11 + invM11*iBlockC.dot(blockC);
    auto temp = -rotateScaleM.rightMultiply(iBlockC);
    VtxCoeffs(vtx, 0, 1) = temp[0];
    VtxCoeffs(vtx, 0, 2) = temp[1];
    VtxCoeffs(vtx, 0, 3) = temp[2];

    VtxCoeffs(vtx, 1, 0) = -iBlockC[0];
    VtxCoeffs(vtx, 1, 1) = invM2D[0] * rotateScaleM(0, 0) + invM2D[1] * rotateScaleM(0, 1) + invM2D[2] * rotateScaleM(0, 2);
    VtxCoeffs(vtx, 1, 2) = invM2D[0] * rotateScaleM(1, 0) + invM2D[1] * rotateScaleM(1, 1) + invM2D[2] * rotateScaleM(1, 2);
    VtxCoeffs(vtx, 1, 3) = invM2D[0] * rotateScaleM(2, 0) + invM2D[1] * rotateScaleM(2, 1) + invM2D[2] * rotateScaleM(2, 2);
   
    VtxCoeffs(vtx, 2, 0) = -iBlockC[1];
    VtxCoeffs(vtx, 2, 1) = invM2D[1] * rotateScaleM(0, 0) + invM2D[3] * rotateScaleM(0, 1) + invM2D[4] * rotateScaleM(0, 2);
    VtxCoeffs(vtx, 2, 2) = invM2D[1] * rotateScaleM(1, 0) + invM2D[3] * rotateScaleM(1, 1) + invM2D[4] * rotateScaleM(1, 2);
    VtxCoeffs(vtx, 2, 3) = invM2D[1] * rotateScaleM(2, 0) + invM2D[3] * rotateScaleM(2, 1) + invM2D[4] * rotateScaleM(2, 2);

    //Calculate a smmooth flag for if we are near a material edge (from MPAS)
    double pMassGradNorm = vtxMatrices(vtx, 1) * vtxMatrices(vtx, 1) +  vtxMatrices(vtx, 2) * vtxMatrices(vtx, 2) +
                           vtxMatrices(vtx, 3) * vtxMatrices(vtx, 3);
    pMassGradNorm = (sqrt(pMassGradNorm) / vtx_area_sqrt) / Kokkos::max(vtxMatrices(vtx, 0), 1e-4);

    double ramp = 2.2 - 10.0 * pMassGradNorm;
    ramp = ramp < 0.0 ? 0.0 : ramp;
    ramp = ramp > 1.0 ? 1.0 : ramp;

    double massRamp = 4.0 * (vtxMatrices(vtx, 0) - 1.0);
    massRamp = massRamp < 0.0 ? 0.0 : massRamp;
    massRamp = massRamp > 1.0 ? 1.0 : massRamp;

    ramp *= massRamp;
    nearAnEdge_l(vtx) = (interiorVertex(vtx) == 1) ? ramp : 0;
  });
  this->precomputedVtxCoeffs_new = VtxCoeffs;
  this->nearAnEdge = nearAnEdge_l;
}

template <MeshFieldIndex meshFieldIndex>
void MPMesh::assemblyVtx1(){
  Kokkos::Timer timer;

  int self, numProcsTot;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  MPI_Comm_size(comm, &numProcsTot);

  auto VtxCoeffs_new=this->precomputedVtxCoeffs_new;

  //Mesh Information
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  int numVtx = p_mesh->getNumVertices();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  int numVertices = p_mesh->getNumVertices();

  //Mesh Field
  constexpr MaterialPointSlice mpfIndex = meshFieldIndexToMPSlice<meshFieldIndex>;
  const int numEntries = mpSliceToNumEntries<mpfIndex>();
  p_mesh->fillMeshField<meshFieldIndex>(numVtx, numEntries, 0.0);
  auto meshField = p_mesh->getMeshField<meshFieldIndex>();

  //Material Points
  auto mpData = p_MPs->getData<mpfIndex>();
  auto weight = p_MPs->getData<MPF_Basis_Vals>();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();

  //Earth Radius
  double radius = 1.0;
  if(p_mesh->getGeomType() == geom_spherical_surf)
    radius=p_mesh->getSphereRadius();

  //Reconstruct
  auto reconstruct = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1;
        double w_vtx=weight(mp,i); 
        double CoordDiffs[vec4d_nEntries] = {1, (-vtxCoords(vID,0) + mpPositions(mp,0))/radius,
                                                (-vtxCoords(vID,1) + mpPositions(mp,1))/radius,
                                                (-vtxCoords(vID,2) + mpPositions(mp,2))/radius};

        auto factor = w_vtx*(VtxCoeffs_new(vID,0, 0) + VtxCoeffs_new(vID,0, 1)*CoordDiffs[1] +
                                                       VtxCoeffs_new(vID,0, 2)*CoordDiffs[2] +
                                                       VtxCoeffs_new(vID,0, 3)*CoordDiffs[3]);

        for (int k=0; k<numEntries; k++){
          auto val = factor*mpData(mp,k);
          Kokkos::atomic_add(&meshField(vID,k), val);
        }
      }
    }
  };
  p_MPs->parallel_for(reconstruct, "reconstruct");
  Kokkos::fence();
  pumipic::RecordTime("Assemble Field per process" + std::to_string(self), timer.seconds());

  timer.reset();
  if(numProcsTot>1){ 
    communicate_and_take_halo_contributions(meshField, numVertices, numEntries, 0, 0);
    communicate_and_take_halo_contributions(meshField, numVertices, numEntries, 1, 1);
  }
  pumipic::RecordTime("Communicate Field Values" + std::to_string(self), timer.seconds());
}

//Start Communication routine
void MPMesh::startCommunication(){
  std::cout<<__FUNCTION__<<std::endl;
  Kokkos::Timer timer;
  int self, numProcsTot;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  MPI_Comm_size(comm, &numProcsTot); 

  //The routine should work for elements too, although currently the communication 
  //is done for vertices. For elements, the follwoing three variables should correspond 
  //to elements.
  auto entOwners = p_mesh->getVtx2Process();
  auto ent2global = p_mesh->getVtxGlobal();
  int numEntities = p_mesh->getNumVertices();

  //Loop over elements and find no of owners and halos
  Kokkos::View<int> owner_count("owner_count");
  Kokkos::View<int> halo_count("halo_count");
  Kokkos::deep_copy(owner_count, 0);
  Kokkos::deep_copy(halo_count, 0);
  Kokkos::parallel_for("countOwnerHalo", numEntities, KOKKOS_LAMBDA(const int elm){
    if (entOwners(elm)==self)
      Kokkos::atomic_add(&owner_count(), 1);
    else
      Kokkos::atomic_add(&halo_count(), 1);
  });

  Kokkos::deep_copy(numOwnersTot, owner_count);
  Kokkos::deep_copy(numHalosTot, halo_count);
  assert(numHalosTot+numOwnersTot == numEntities);
  int num_ints_per_copy = 2;

  //#Halo Cells/proc which are owners on other process
  numOwnersOnOtherProcs.resize(numProcsTot);

  //#OwnerCells/proc which are halos on other proces
  numHalosOnOtherProcs.resize(numProcsTot); 

  //For every halo cell find the owning process and the local Id in that process
  haloOwnerProcs.reserve(numHalosTot);
  haloOwnerLocalIDs.resize(numProcsTot);

  //Copy owning processes and globalIds to CPU
  auto entOwners_host = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace::memory_space(),
                        entOwners);
  auto ent2global_host = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace::memory_space(),
                         ent2global);

  //Do Map of Global To Local ID
  //Check unordered vs ordered map; which faster?
  std::map<int, int> global2local;
  //std::unordered_map<int, int> global2local;
  for (int iEnt = 0; iEnt < numEntities; iEnt++) {
    int globalID = ent2global_host(iEnt);
    global2local[globalID] = iEnt;
  }

  //Loop over all halo Entities and find the owning process
  for (auto iEnt=numOwnersTot; iEnt<numOwnersTot+numHalosTot; iEnt++){
    auto ownerProc = entOwners_host[iEnt];
    assert(entOwners_host(iEnt) != self);
    numOwnersOnOtherProcs[ownerProc] = numOwnersOnOtherProcs[ownerProc]+1;
    haloOwnerProcs.push_back(ownerProc);
  }

  MPI_Alltoall(numOwnersOnOtherProcs.data(), 1, MPI_INT, numHalosOnOtherProcs.data(), 1, MPI_INT, comm);

  // Halo Entity's Global & Local Id To Owning Process
  std::vector<std::vector<int>> sendBufs(numProcsTot);
  for (int proc = 0; proc < numProcsTot; proc++)
    sendBufs[proc].reserve(num_ints_per_copy*numOwnersOnOtherProcs[proc]);

  for (int iEnt=numOwnersTot; iEnt<numOwnersTot+numHalosTot; iEnt++) {
    auto ownerProc = entOwners_host(iEnt);
    assert(ownerProc != self);
    sendBufs[ownerProc].push_back(ent2global_host(iEnt));
    sendBufs[ownerProc].push_back(iEnt);
  }

  //Requests  
  std::vector<MPI_Request> requests;
  requests.reserve(2*numProcsTot);

  //Receive Calls
  std::vector<std::vector<int>> recvBufs(numProcsTot);
  for (int proc = 0; proc < numProcsTot; proc++) {
    if (numHalosOnOtherProcs[proc] > 0) {
      recvBufs[proc].resize(num_ints_per_copy*numHalosOnOtherProcs[proc]);
      MPI_Request req;
      MPI_Irecv(recvBufs[proc].data(), num_ints_per_copy*numHalosOnOtherProcs[proc], MPI_INT, proc, MPI_ANY_TAG, comm, &req);
      requests.push_back(req);
    }
  }
  //Send Calls
  for (int proc=0; proc<numProcsTot; proc++) {
    auto& buf=sendBufs[proc];
    if(buf.empty()) continue;
    MPI_Request req;
    MPI_Isend(buf.data(), buf.size(), MPI_INT, proc, 0, comm, &req);
    requests.push_back(req);
  }
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
  requests.clear();

  //Now the owner process needs to look at these globalIDs convert them to localIds and send it back
  //recvBufs[p] contains global IDs of elements that halo rank p needs
  //numHalosOnOtherProcs[p] tells how many to expect from proc p
  ownerOwnerLocalIDs.resize(numProcsTot);
  ownerHaloLocalIDs.resize(numProcsTot);

  for (int proc = 0; proc < numProcsTot; proc++) {
    if (numHalosOnOtherProcs[proc] > 0) {
      ownerOwnerLocalIDs[proc].resize(numHalosOnOtherProcs[proc]);
      ownerHaloLocalIDs[proc].resize(numHalosOnOtherProcs[proc]);
      for (int i = 0; i < numHalosOnOtherProcs[proc]; i++) {
        int globalID = recvBufs[proc][i*num_ints_per_copy];
        ownerOwnerLocalIDs[proc][i] = global2local[globalID];
        ownerHaloLocalIDs[proc][i]  = recvBufs[proc][i*num_ints_per_copy+1];
      }
    }
  }

  // On the halo side, need to receive the localIds of owning Process
  for (int proc = 0; proc < numProcsTot; proc++) {
    if (numOwnersOnOtherProcs[proc] > 0) { // these are cells whose owners are in other processes
      haloOwnerLocalIDs[proc].resize(numOwnersOnOtherProcs[proc]);
      MPI_Request req;
      MPI_Irecv(haloOwnerLocalIDs[proc].data(), haloOwnerLocalIDs[proc].size(), MPI_INT, proc, MPI_ANY_TAG, comm, &req);
      requests.push_back(req);
    }
  }

  //Sends back localID of the owned cells so that HaloToOwner can be done for halo processes
  for (int proc = 0; proc < numProcsTot; proc++) {
    if (numHalosOnOtherProcs[proc]>0) {
      MPI_Request req;
      MPI_Isend(ownerOwnerLocalIDs[proc].data(), ownerOwnerLocalIDs[proc].size(), MPI_INT, proc, 1, comm, &req);
      requests.push_back(req);
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

  pumipic::RecordTime("Start Communication" + std::to_string(self), timer.seconds());
}

void MPMesh::communicate_and_take_halo_contributions(const Kokkos::View<double**>& meshField, int nEntities, int numEntries, int mode, int op){
  int self;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);

  Kokkos::Timer timer;
  auto reconVals_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), meshField);
  Kokkos::fence();
  std::vector<std::vector<double>> fieldData(nEntities, std::vector<double>(numEntries, 0.0));
  for (int i = 0; i < nEntities; ++i) {
    for (int j = 0; j < numEntries; ++j) {
      fieldData[i][j] = reconVals_host(i, j);
    }
  }
  pumipic::RecordTime("Communication-GPU to CPU-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());

  timer.reset();
  std::vector<std::vector<int>>    recvIDVec;
  std::vector<std::vector<double>> recvDataVec;
  communicateFields(fieldData, nEntities, numEntries, mode, recvIDVec, recvDataVec);
  pumipic::RecordTime("Communication-InterProcess-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());

  timer.reset();
  int numProcsTot =  recvIDVec.size();
  //Flatten IDs 
  int totalSize = 0;
  std::vector<int> offsets(numProcsTot, 0); 
  for(int i=0; i<numProcsTot; i++) {
    offsets[i] = totalSize;
    totalSize += recvIDVec[i].size();
  }
  std::vector<int> flatIDVec(totalSize, 0);
  for(int i=0; i<numProcsTot; i++) {
    std::copy(recvIDVec[i].begin(), recvIDVec[i].end(), flatIDVec.begin() + offsets[i]);
  }
  Kokkos::View<int*> recvIDGPU("recvIDGPU", totalSize);
  auto hostView = Kokkos::View<int*, Kokkos::HostSpace>("recvIDCPU", totalSize);
  std::copy(flatIDVec.begin(), flatIDVec.end(), hostView.data());
  Kokkos::deep_copy(recvIDGPU, hostView);

  //Flatten Data
  int totalSize_data=0;
  std::vector<int> offsets_data(numProcsTot, 0);
  for(int i=0; i<numProcsTot; i++) {
    offsets_data[i] = totalSize_data;
    totalSize_data += recvDataVec[i].size();
  }
  std::vector<double> flatDataVec(totalSize_data, 0);
  for(int i=0; i<numProcsTot; i++) {
    std::copy(recvDataVec[i].begin(), recvDataVec[i].end(), flatDataVec.begin() + offsets_data[i]);
  }
  Kokkos::View<double*> recvDataGPU("recvDataGPU", totalSize_data);
  auto hostView_data= Kokkos::View<double*, Kokkos::HostSpace>("recvDataCPU", totalSize_data);
  std::copy(flatDataVec.begin(), flatDataVec.end(), hostView_data.data()); 
  Kokkos::deep_copy(recvDataGPU, hostView_data);
  Kokkos::fence();
  //Assertions
  assert(totalSize_data == totalSize*numEntries);
  for (int i=0; i<numProcsTot; i++){
    assert(recvDataVec[i].size() == recvIDVec[i].size() * numEntries);
  }
  pumipic::RecordTime("Communication-CPU to GPU-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());
  
  //Take contributions from other procs
  timer.reset();
  Kokkos::parallel_for("halo contribution", recvIDGPU.size(), KOKKOS_LAMBDA(const int i){
    int vertex = recvIDGPU(i);
    for(int k=0; k<numEntries; k++){
      if(op==0) Kokkos::atomic_add(&meshField(vertex,k), recvDataGPU(i*numEntries+k));
      if(op==1) meshField(vertex, k) = recvDataGPU(i * numEntries + k);
    }
  });
  Kokkos::fence();
  pumipic::RecordTime("Communication-GPU reduction-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());
}

void MPMesh::communicateFields(const std::vector<std::vector<double>>& fieldData, const int numEntities, const int numEntries, int mode, 
                               std::vector<std::vector<int>>& recvIDVec,  std::vector<std::vector<double>>& recvDataVec){
  int self, numProcsTot;
  MPI_Comm comm = p_MPs->getMPIComm();
  MPI_Comm_rank(comm, &self);
  MPI_Comm_size(comm, &numProcsTot);

  assert(numEntities == numOwnersTot + numHalosTot);

  std::vector<std::vector<double>> sendDataVec(numProcsTot);

  recvIDVec.resize(numProcsTot);
  recvDataVec.resize(numProcsTot);

  for(int i = 0; i < numProcsTot; i++){
    if(i==self) continue;

    int numToSend = 0, numToRecv = 0; 
    if(mode == 0) {
      //gather (halos send to owners)
      numToSend = numOwnersOnOtherProcs[i];
      numToRecv = numHalosOnOtherProcs[i];
    }
    else{ 
      //scatter (owners send to halos)
      numToSend = numHalosOnOtherProcs[i];
      numToRecv = numOwnersOnOtherProcs[i];
    }

    if(numToSend > 0){
      sendDataVec[i].reserve(numToSend*numEntries);
    }
    if(numToRecv > 0){
      recvDataVec[i].resize(numToRecv*numEntries);
      recvIDVec[i].resize(numToRecv);
    }
  }

  if(mode == 0){
    // Halos sends to owners
    for (int iEnt = 0; iEnt < numHalosTot; iEnt++){
      auto ownerProc = haloOwnerProcs[iEnt];
      for (int iDouble = 0; iDouble < numEntries; iDouble++)
        sendDataVec[ownerProc].push_back(fieldData[numOwnersTot+iEnt][iDouble]);
    }
  }
  else if(mode == 1){
    // Owner sends to halos
    for (size_t iProc=0; iProc<ownerOwnerLocalIDs.size(); iProc++) {
      for (auto& ownerID : ownerOwnerLocalIDs[iProc]) {
        for (int iDouble = 0; iDouble < numEntries; iDouble++)
          sendDataVec[iProc].push_back(fieldData[ownerID][iDouble]);
      }
    }
  }

  std::vector<MPI_Request> requests;
  requests.reserve(4*numProcsTot);
  for(int proc = 0; proc < numProcsTot; proc++){
    if(proc == self) continue;  
    if(mode == 0 && numHalosOnOtherProcs[proc]){
      assert(recvIDVec[proc].size() == (size_t)numHalosOnOtherProcs[proc]);
      assert(recvDataVec[proc].size() == recvIDVec[proc].size() * (size_t)numEntries);
      MPI_Request req3, req4;
      MPI_Irecv(recvIDVec[proc].data(), recvIDVec[proc].size(), MPI_INT, proc, 1, comm, &req3);
      MPI_Irecv(recvDataVec[proc].data(), recvDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req4);
      requests.push_back(req3);
      requests.push_back(req4);
    }
    if(mode == 0 && numOwnersOnOtherProcs[proc]) {
      assert(haloOwnerLocalIDs[proc].size() == (size_t)numOwnersOnOtherProcs[proc]);
      assert(sendDataVec[proc].size() == haloOwnerLocalIDs[proc].size() * (size_t)numEntries);
      MPI_Request req1, req2;
      MPI_Isend(haloOwnerLocalIDs[proc].data(), haloOwnerLocalIDs[proc].size(), MPI_INT, proc, 1, comm, &req1);
      MPI_Isend(sendDataVec[proc].data(), sendDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req2);
      requests.push_back(req1);
      requests.push_back(req2);
    }

    if(mode == 1 && numOwnersOnOtherProcs[proc]){
      MPI_Request req3, req4;
      MPI_Irecv(recvIDVec[proc].data(), recvIDVec[proc].size(), MPI_INT, proc, 1, comm, &req3);
      MPI_Irecv(recvDataVec[proc].data(), recvDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req4);
      requests.push_back(req3);
      requests.push_back(req4);
    }
    if(mode == 1 && numHalosOnOtherProcs[proc]) {
      MPI_Request req1, req2;
      MPI_Isend(ownerHaloLocalIDs[proc].data(), ownerHaloLocalIDs[proc].size(), MPI_INT, proc, 1, comm, &req1);
      MPI_Isend(sendDataVec[proc].data(), sendDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req2);
      requests.push_back(req1);
      requests.push_back(req2);
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

template <MeshFieldIndex meshFieldIndex>
void MPMesh::assembly(int order, MeshFieldType type, bool basisWeightFlag, bool massWeightFlag){
  if(basisWeightFlag || massWeightFlag) {
    std::cerr << "WARNING: basis and mass weight flags ignored\n";
  }

  if (order == 0 && type == MeshFType_VtxBased)
    assemblyVtx0<meshFieldIndex>();
  else if (order == 0 && type == MeshFType_ElmBased)
    assemblyElm0<meshFieldIndex>();
  else if (order == 1 && type == MeshFType_VtxBased)
    assemblyVtx1<meshFieldIndex>();
  else{
    std::cerr << "Error: Assembly order is not supported\n";
    exit(1);
  }
}

// (HDT) weighted assembly of scalar field
template <MaterialPointSlice index>
DoubleView MPMesh::wtScaAssembly(){
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  int numVtxs = p_mesh->getNumVertices(); // total number of vertices of the mesh
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();

  DoubleView vField("wtScaField", numVtxs); // Kokkos array of double type, size = numVtxs

  auto mpData = p_MPs->getData<index>();

  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if (mask) {
      /* get the coordinates of all the vertices of elm */
      int nElmVtxs = elm2VtxConn(elm,0);      // number of vertices bounding the element
      Vec2d eVtxCoords[maxVtxsPerElm + 1];
      for (int i = 1; i <= nElmVtxs; i++) {
        // elm2VtxConn(elm,i) is the vertex ID (1-based index) of vertex #i of elm
        eVtxCoords[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
        eVtxCoords[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
      }
      // last component of eVtxCoords stores the firs vertex (to avoid if-condition in the Wachspress computation)
      eVtxCoords[nElmVtxs][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
      eVtxCoords[nElmVtxs][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
      
      /* compute the values of basis functions at mp position */
      double basisByArea[maxElmsPerVtx];
      Vec2d mpCoord(mpPositions(mp,0), mpPositions(mp,1));
      getBasisByAreaGblForm(mpCoord, nElmVtxs, eVtxCoords, basisByArea);

      /* get the mp's property that is assembled to vertices */
      double assemVal = mpData(mp, 0); // ??? for scalar mp data, is index 0 always?

      /* accumulate the mp's property to vertices */
      for (int i = 0; i < nElmVtxs; i++) {
        int vID = elm2VtxConn(elm,i+1)-1;
        Kokkos::atomic_add(&vField(vID), assemVal * basisByArea[i]);
      }
    }
  };
  p_MPs->parallel_for(assemble, "weightedScalarAssembly");
  return vField;
} // wtScaAssembly

// (HDT) weighted assembly of vector2 field (not weighted by mass/volume)
template <MaterialPointSlice index>
Vec2dView MPMesh::wtVec2Assembly(){
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  int numVtxs = p_mesh->getNumVertices(); // total number of vertices of the mesh
  auto elm2VtxConn = p_mesh->getElm2VtxConn();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();

  Vec2dView vField("wtVec2Field", numVtxs); // Kokkos array of Vec2d type, size = numVtxs

  auto mpData = p_MPs->getData<index>();

  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if (mask) {
      /* collect the coordinates of all the vertices of elm */
      int nElmVtxs = elm2VtxConn(elm,0);      // number of vertices bounding the element
      Vec2d eVtxCoords[maxVtxsPerElm + 1];
      for (int i = 1; i <= nElmVtxs; i++) {
        // elm2VtxConn(elm,i) is the vertex ID (1-based index) of vertex #i of elm
        eVtxCoords[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);    
        eVtxCoords[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
      }
      // last component of eVtxCoords stores the firs vertex (to avoid if-condition in the Wachspress computation)
      eVtxCoords[nElmVtxs][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
      eVtxCoords[nElmVtxs][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);

      /* compute the values of basis functions at mp position */
      double basisByArea[maxElmsPerVtx];
      Vec2d mpCoord(mpPositions(mp,0), mpPositions(mp,1));
      getBasisByAreaGblForm(mpCoord, nElmVtxs, eVtxCoords, basisByArea);

      /* get the mp's volume */
      double mpVolume = 1.0; // TODO: change to mp's volume here

      /* get the mp's property to be assembled */
      Vec2d assemVal;
      assemVal[0] =  mpData(mp, 0) * mpVolume;
      assemVal[1] =  mpData(mp, 1) * mpVolume;

      /* accumulate the mp's constructed quantities to the cell vertices */
      for (int i = 0; i < nElmVtxs; i++) {
        int vID = elm2VtxConn(elm,i+1)-1;
        Kokkos::atomic_add(&(vField(vID)[0]), assemVal[0] * basisByArea[i]);
        Kokkos::atomic_add(&(vField(vID)[1]), assemVal[1] * basisByArea[i]);
      }
    }
  };
  p_MPs->parallel_for(assemble, "weightedVec2dAssembly");
  return vField;
} // wtVec2Assembly

template<MeshFieldIndex meshFieldIndex>
void MPMesh::setReconstructSlice(int order, MeshFieldType type) {
  auto function = [=](){ assembly<meshFieldIndex>(order, type, false, false); };
  const auto [iter, success] = reconstructSlice.insert({meshFieldIndex, function});
  if (!success){
    std::cerr << "Error: Slice is already being reconstructed\n";
    exit(1);
  }
}

} //end namespace polyMPO
#endif
