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
    //for elm in elementsInMesh { //pseudo code - the 'parallel_for' handles this
    //  for mp in materialPointsInElm { //pseudo code (cont.)
          if(mask) { //if material point is 'active'/'enabled'
            int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
            for(int i=0; i<nVtxE; i++){
              int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
              double distance = mpPositions(mp,0) + mpPositions(mp,1) + mpPositions(mp,2);
              Kokkos::atomic_add(&vField(vID),distance);
            }
          }
    //  }
    //}
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
  
  std::cout<<__FUNCTION__<<std::endl;	
  Kokkos::Timer timer;
  constexpr MaterialPointSlice mpfIndex = meshFieldIndexToMPSlice<meshFieldIndex>;
  auto mpData = p_MPs->getData<mpfIndex>();
  const int numEntries = mpSliceToNumEntries<mpfIndex>();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();

  int numElms = p_mesh->getNumElements();
  p_mesh->fillMeshField<meshFieldIndex>(numElms, numEntries, 0.0);
  auto meshField = p_mesh->getMeshField<meshFieldIndex>();

  Kokkos::View<int*> mpsPerElm("mpsPerElm", numElms);
  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      Kokkos::atomic_add(&mpsPerElm(elm),1);
      for(int j=0;j<numEntries;j++){
        Kokkos::atomic_add(&meshField(elm,j), mpData(mp,0));
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

template <MeshFieldIndex meshFieldIndex>
void MPMesh::assemblyVtx1() {

  std::cout<<__FUNCTION__<<std::endl;
  static int count=0;  
  double epsilon=1e-13;  
  //Mesh Inforamtion
  auto elm2VtxConn = p_mesh->getElm2VtxConn();  
  int numVtx = p_mesh->getNumVertices();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();

  //Mesh Field
  constexpr MaterialPointSlice mpfIndex = meshFieldIndexToMPSlice<meshFieldIndex>;
  const int numEntries = mpSliceToNumEntries<mpfIndex>();
  p_mesh->fillMeshField<meshFieldIndex>(numVtx, numEntries, 0.0);
  auto meshField = p_mesh->getMeshField<meshFieldIndex>();

  //Material Points
  auto mpData = p_MPs->getData<mpfIndex>();
  auto weight = p_MPs->getData<MPF_Basis_Vals>();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();

  //Matrix for each vertex
  Kokkos::View<double*[4][4]> VtxMatrices("VtxMatrices", p_mesh->getNumVertices());

  //Reconstructed values
  Kokkos::View<double*> reconVals("meshField", p_mesh->getNumVertices());
  //Kokkos::deep_copy(reconVals, 0.0);
  Kokkos::View<int*>mps_around_vertex("MPs_around_vertex", p_mesh->getNumVertices());
  Kokkos::View<int*>ill_cond_vertex("ill_cond", p_mesh->getNumVertices());
  

  //Earth Radius
  double radius = p_mesh->getSphereRadius();

  //Assemble matrix for each vertex
  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
        double w_vtx=weight(mp,i);
        
        double CoordDiffs[4] = {1, (vtxCoords(vID,0) - mpPositions(mp,0))/radius, (vtxCoords(vID,1) - mpPositions(mp,1))/radius, 
          (vtxCoords(vID,2) - mpPositions(mp,2))/radius};		            //add to the matrix
        for (int k=0; k<4; k++)
          for(int l=0;l<4;l++)
            Kokkos::atomic_add(&VtxMatrices(vID,k,l), CoordDiffs[k] * CoordDiffs[l] * w_vtx);
	
	Kokkos::atomic_add(&mps_around_vertex[vID],1);
      }
    }
  };
  p_MPs->parallel_for(assemble, "assembly");
  
  //Solve Ax=b for each vertex  
  printf("=====Count %d=========\n", count);
  Kokkos::View<double*[4]> VtxCoeffs("VtxMatrices", p_mesh->getNumVertices());
  Kokkos::parallel_for("solving Ax=b", numVtx, KOKKOS_LAMBDA(const int vtx){
    Vec4d v0 = {VtxMatrices(vtx,0,0), VtxMatrices(vtx,0,1), VtxMatrices(vtx,0,2), VtxMatrices(vtx,0,3)};
    Vec4d v1 = {VtxMatrices(vtx,1,0), VtxMatrices(vtx,1,1), VtxMatrices(vtx,1,2), VtxMatrices(vtx,1,3)};
    Vec4d v2 = {VtxMatrices(vtx,2,0), VtxMatrices(vtx,2,1), VtxMatrices(vtx,2,2), VtxMatrices(vtx,2,3)};
    Vec4d v3 = {VtxMatrices(vtx,3,0), VtxMatrices(vtx,3,1), VtxMatrices(vtx,3,2), VtxMatrices(vtx,3,3)};
    Matrix A = {v0,v1,v2,v3};
  
    Vec4d v000={0, 0, 0, 0}; 
    Matrix Q={v000, v000, v000, v000};
    Matrix R={v000, v000, v000, v000};
    QR_decomp(A, Q, R);
     
    /* 
    //For checking QR on a  non illconditioned matrix
    Vec4d v00={1, 0.2, 0.3, 0.4};
    Vec4d v01={0.1, 1, 0.1, 0.1};
    Vec4d v02={0.1, 0.1, 1, 0.1};
    Vec4d v03={0.1, 0.1, 0.1, 1};
    Matrix A_new={v00, v01, v02, v03}; 
    QR_decomp(A,Q,R);
    if(vtx==2050){
      for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
          printf("%.15e  ", Q(i,j));    
        }
	printf ("\n");
      }
     
      for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
          printf("%.15e  ", R(i,j));    
        }
        printf ("\n");
      }
    }
    */
    if(R(0,0)<epsilon || R(1,1)<epsilon || R(2,2)<epsilon || R(3,3)<epsilon) ill_cond_vertex(vtx)=1;


    double coeff[4]={0.0, 0.0, 0.0, 0.0};
    CholeskySolve(A, coeff);
    for (int i=0; i<4; i++) VtxCoeffs(vtx,i)=coeff[i];
    
    if(vtx ==2050){
      printf("Matrix %.15e %.15e %.15e %.15e\n",VtxMatrices(vtx,0,0),VtxMatrices(vtx,0,1),VtxMatrices(vtx,0,2),VtxMatrices(vtx,0,3));
      printf("Matrix %.15e %.15e %.15e \n", VtxMatrices(vtx,1,1),VtxMatrices(vtx,1,2),VtxMatrices(vtx,1,3));
      printf("Matrix %.15e %.15e \n", VtxMatrices(vtx,2,2),VtxMatrices(vtx,2,3));
      printf("Matrix %.15e \n", VtxMatrices(vtx,3,3));
      printf("Vtx %d coeff %.15e %.15e %.15e %.15e \n", vtx, coeff[0], coeff[1], coeff[2], coeff[3]);
    }
    /*
    if(mps_around_vertex(vtx)<4)
      printf("Vtx %d Ill condition number %d and no of surrounding MPs %d \n", vtx, ill_cond_vertex(vtx), mps_around_vertex(vtx));
    if(mps_around_vertex(vtx)>=4 && ill_cond_vertex(vtx)==1)
      printf("Vertex %d Anomaly \n", vtx);
    */

  });
 
  //Reconstruct
  auto reconstruct = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1;
        double w_vtx=weight(mp,i); 
        double CoordDiffs[4] = {1, (vtxCoords(vID,0) - mpPositions(mp,0))/radius, (vtxCoords(vID,1) - mpPositions(mp,1))/radius, 
          (vtxCoords(vID,2) - mpPositions(mp,2))/radius};
	auto val = w_vtx*(VtxCoeffs(vID,0) + VtxCoeffs(vID,1)*CoordDiffs[1] + VtxCoeffs(vID,2)*CoordDiffs[2] + 
	  VtxCoeffs(vID,3)*CoordDiffs[3])*mpData(mp,0) ;
	Kokkos::atomic_add(&reconVals(vID), val);
        if(vID==2050){
           printf("Vtx %d coeff val %.15e value %.15e\n", vID, val, reconVals(vID));
        }
      }
    }
  };
  p_MPs->parallel_for(reconstruct, "reconstruct");
  
  Kokkos::parallel_for("averaging", numVtx, KOKKOS_LAMBDA(const int vtx){
    if(mps_around_vertex(vtx)<4) reconVals(vtx)=-1;
    if(ill_cond_vertex(vtx)==1)  reconVals(vtx)=-2; 
    meshField(vtx, 0) = reconVals(vtx);
  });
  count ++;
  
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

        /* get the mp's property that is assebled to vertices */
        double assValue = mpData(mp, 0); // ??? for scalar mp data, is index 0 always?

        /* accumulate the mp's property to vertices */
        for (int i = 0; i < nElmVtxs; i++) {
          int vID = elm2VtxConn(elm,i+1)-1;
          Kokkos::atomic_add(&vField(vID), assValue * basisByArea[i]);
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
        Vec2d assValue;
        assValue[0] =  mpData(mp, 0) * mpVolume;
        assValue[1] =  mpData(mp, 1) * mpVolume;

        /* accumulate the mp's constructed quantities to the cell vertices */
        for (int i = 0; i < nElmVtxs; i++) {
          int vID = elm2VtxConn(elm,i+1)-1;
          Kokkos::atomic_add(&(vField(vID)[0]), assValue[0] * basisByArea[i]);
          Kokkos::atomic_add(&(vField(vID)[1]), assValue[1] * basisByArea[i]);
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
