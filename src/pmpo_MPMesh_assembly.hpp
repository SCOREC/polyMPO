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


void MPMesh::resetPreComputeFlag(){
  isPreComputed = false;
}

void MPMesh::computeMatricesAndSolve(){
  Kokkos::Timer timer;
  //Mesh Information
  auto elm2VtxConn = p_mesh->getElm2VtxConn();  
  int numVtx = p_mesh->getNumVertices();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();

  //Dual Element Area for Regularization
  auto dual_triangle_area=p_mesh->getMeshField<MeshF_DualTriangleArea>();

  //Material Points
  auto weight = p_MPs->getData<MPF_Basis_Vals>();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();

  //Matrix for each vertex
  Kokkos::View<double*[vec4d_nEntries][vec4d_nEntries]> VtxMatrices("VtxMatrices", p_mesh->getNumVertices());

  //Earth Radius
  double radius = 1.0;
  if(p_mesh->getGeomType() == geom_spherical_surf)
    radius=p_mesh->getSphereRadius();

  bool scaling=true;
  int reg_method = 2;
  
  //Assemble matrix for each vertex
  auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
        double w_vtx=weight(mp,i);
        double mScale=1;
        if(scaling)
          mScale=sqrt(dual_triangle_area(vID,0))/radius;
        double CoordDiffs[vec4d_nEntries] = {1, (vtxCoords(vID,0) - mpPositions(mp,0))/radius, 
                                                (vtxCoords(vID,1) - mpPositions(mp,1))/radius, 
                                                (vtxCoords(vID,2) - mpPositions(mp,2))/radius};		 
        //All entries except first row and column
        for (int k=1; k<vec4d_nEntries; k++)
          for (int l=1; l<vec4d_nEntries; l++)
            Kokkos::atomic_add(&VtxMatrices(vID,k,l), CoordDiffs[k] * CoordDiffs[l] * w_vtx);
	//First entry
        Kokkos::atomic_add(&VtxMatrices(vID,0,0), CoordDiffs[0] * CoordDiffs[0] * w_vtx*mScale*mScale); 
	//First row and column except the first entry
	for (int k=1; k<vec4d_nEntries; k++){
          Kokkos::atomic_add(&VtxMatrices(vID,0,k), CoordDiffs[0] * CoordDiffs[k] * w_vtx*mScale);
	  Kokkos::atomic_add(&VtxMatrices(vID,k,0), CoordDiffs[k] * CoordDiffs[0] * w_vtx*mScale);
	}
      }
    }
  };
  p_MPs->parallel_for(assemble, "assembly");
  
  //Assemble matrix for each vertex and apply regularization
  Kokkos::View<double*[vec4d_nEntries]> VtxCoeffs("VtxCoeffs", p_mesh->getNumVertices());

  Kokkos::parallel_for("solving Ax=b", numVtx, KOKKOS_LAMBDA(const int vtx){
    Vec4d v0 = {VtxMatrices(vtx,0,0), VtxMatrices(vtx,0,1), VtxMatrices(vtx,0,2), VtxMatrices(vtx,0,3)};
    Vec4d v1 = {VtxMatrices(vtx,1,0), VtxMatrices(vtx,1,1), VtxMatrices(vtx,1,2), VtxMatrices(vtx,1,3)};
    Vec4d v2 = {VtxMatrices(vtx,2,0), VtxMatrices(vtx,2,1), VtxMatrices(vtx,2,2), VtxMatrices(vtx,2,3)};
    Vec4d v3 = {VtxMatrices(vtx,3,0), VtxMatrices(vtx,3,1), VtxMatrices(vtx,3,2), VtxMatrices(vtx,3,3)};
    //Define the matrices
    Matrix4d A = {v0,v1,v2,v3};
    Matrix4d A_regularized = {v0, v1, v2, v3};
    //Regularization
    switch(reg_method){
      case 0:{
        break;
      }   
      case 1:{
        double A_trace = A.trace();
        A_regularized.addToDiag(A_trace*1e-8);
        break;
      }
      case 2:{
        double regParam=sqrt(EPSILON)*(VtxMatrices(vtx,0,0)+VtxMatrices(vtx,1,1)+
			               VtxMatrices(vtx,2,2)+VtxMatrices(vtx,3,3));
        A_regularized.addToDiag(regParam);
        break;
      }
      default:{
        printf("Invalid regularization method \n");
        break;	
      }
    }
    //Solve Ax=b 
    double coeff[vec4d_nEntries]={0.0, 0.0, 0.0, 0.0};
    CholeskySolve4d_UnitRHS(A_regularized, coeff);
    // Undo scaling
    double mScale=1;
    if(scaling)
      mScale=sqrt(dual_triangle_area(vtx,0))/radius;
        
    coeff[0]=coeff[0]*mScale*mScale;
    coeff[1]=coeff[1]*mScale;
    coeff[2]=coeff[2]*mScale;
    coeff[3]=coeff[3]*mScale;
    for (int i=0; i<vec4d_nEntries; i++) 
      VtxCoeffs(vtx,i)=coeff[i];
  });
  this->precomputedVtxCoeffs = VtxCoeffs;
  pumipic::RecordTime("PolyMPO_Calculate_MLS_Coeff", timer.seconds());
}

template <MeshFieldIndex meshFieldIndex>
void MPMesh::assemblyVtx1() {
  Kokkos::Timer timer; 
  //If no reconstruction till now calculate the coeffs
  if (!isPreComputed) {
    computeMatricesAndSolve();
    isPreComputed=true;
  }
  
  auto VtxCoeffs=this->precomputedVtxCoeffs;
  //Mesh Information
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

  //Earth Radius
  double radius = 1.0;
  if(p_mesh->getGeomType() == geom_spherical_surf)
    radius=p_mesh->getSphereRadius();

  //Reconstructed values
  Kokkos::View<double**> reconVals("meshField", p_mesh->getNumVertices(), numEntries);

  //Reconstruct
  auto reconstruct = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1;
        double w_vtx=weight(mp,i); 
        double CoordDiffs[vec4d_nEntries] = {1, (vtxCoords(vID,0) - mpPositions(mp,0))/radius,
                                                (vtxCoords(vID,1) - mpPositions(mp,1))/radius, 
                                                (vtxCoords(vID,2) - mpPositions(mp,2))/radius};

        auto factor = w_vtx*(VtxCoeffs(vID,0) + VtxCoeffs(vID,1)*CoordDiffs[1] + 
                                                VtxCoeffs(vID,2)*CoordDiffs[2] + 
                                                VtxCoeffs(vID,3)*CoordDiffs[3]);
  
        for (int k=0; k<numEntries; k++){
          auto val = factor*mpData(mp,k);
          Kokkos::atomic_add(&reconVals(vID,k), val);
        }
      }
    }
  };
  p_MPs->parallel_for(reconstruct, "reconstruct");

  //Assign as a field 
  Kokkos::parallel_for("assigning", numVtx, KOKKOS_LAMBDA(const int vtx){
    for(int k=0; k<numEntries; k++)
      meshField(vtx, k) = reconVals(vtx,k);
  });
  pumipic::RecordTime("PolyMPO_Reconstruct_Vtx1", timer.seconds());
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
