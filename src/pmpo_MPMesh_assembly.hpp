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

// Linear Reconstruction Method 1, no MPI.
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

//Method 2: Uses subassembly, depends on MPAS for MPI
void MPMesh::subAssemblyCoeffs(int vtxPerElm, int nCells, double* m11, double* m12, double* m13, double* m14, 
                                                          double* m22, double* m23, double* m24, 
                                                          double* m33, double* m34, 
                                                          double* m44){
  
  Kokkos::Timer timer;
  
  MPI_Comm comm = p_MPs->getMPIComm(); 
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);

  //Mesh Information
  auto elm2VtxConn = p_mesh->getElm2VtxConn();  
  int numVtx = p_mesh->getNumVertices();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  auto elm2Process = p_mesh->getElm2Process();
  auto elm2global = p_mesh->getElmGlobal();

  //Dual Element Area for Regularization
  auto dual_triangle_area=p_mesh->getMeshField<MeshF_DualTriangleArea>();
  
  //Material Points
  calcBasis();
  auto weight = p_MPs->getData<MPF_Basis_Vals>();
  auto mpPos = p_MPs->getData<MPF_Cur_Pos_XYZ>();
  auto mpAppID = p_MPs->getData<polyMPO::MPF_MP_APP_ID>();

  //Radius
  double radius = 1.0;
  if(p_mesh->getGeomType() == geom_spherical_surf)
    radius=p_mesh->getSphereRadius();

  Kokkos::View<double**> m11_d("m11", vtxPerElm, nCells);
  Kokkos::View<double**> m12_d("m12", vtxPerElm, nCells);
  Kokkos::View<double**> m13_d("m13", vtxPerElm, nCells);
  Kokkos::View<double**> m14_d("m14", vtxPerElm, nCells);
  Kokkos::View<double**> m22_d("m22", vtxPerElm, nCells);
  Kokkos::View<double**> m23_d("m23", vtxPerElm, nCells);
  Kokkos::View<double**> m24_d("m23", vtxPerElm, nCells);
  Kokkos::View<double**> m33_d("m33", vtxPerElm, nCells);
  Kokkos::View<double**> m34_d("m34", vtxPerElm, nCells);
  Kokkos::View<double**> m44_d("m34", vtxPerElm, nCells);
 
  auto sub_assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask && (elm2Process(elm)==comm_rank)) { //if material point is 'active'/'enabled'
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
          
          double w_vtx=weight(mp,i);
          double mScale=sqrt(dual_triangle_area(vID,0))/radius;
          
          Kokkos::atomic_add(&m11_d(i,elm), w_vtx*mScale*mScale);
          Kokkos::atomic_add(&m12_d(i,elm), w_vtx*mScale*(vtxCoords(vID,0)-mpPos(mp,0))/radius);
          Kokkos::atomic_add(&m13_d(i,elm), w_vtx*mScale*(vtxCoords(vID,1)-mpPos(mp,1))/radius);
          Kokkos::atomic_add(&m14_d(i,elm), w_vtx*mScale*(vtxCoords(vID,2)-mpPos(mp,2))/radius);
          Kokkos::atomic_add(&m22_d(i,elm), w_vtx*mScale*(vtxCoords(vID,0)-mpPos(mp,0))*(vtxCoords(vID,0)-mpPos(mp,0))/(radius*radius));
          Kokkos::atomic_add(&m23_d(i,elm), w_vtx*mScale*(vtxCoords(vID,0)-mpPos(mp,0))*(vtxCoords(vID,1)-mpPos(mp,1))/(radius*radius));
          Kokkos::atomic_add(&m24_d(i,elm), w_vtx*mScale*(vtxCoords(vID,0)-mpPos(mp,0))*(vtxCoords(vID,2)-mpPos(mp,2))/(radius*radius));
          Kokkos::atomic_add(&m33_d(i,elm), w_vtx*mScale*(vtxCoords(vID,1)-mpPos(mp,1))*(vtxCoords(vID,1)-mpPos(mp,1))/(radius*radius));
          Kokkos::atomic_add(&m34_d(i,elm), w_vtx*mScale*(vtxCoords(vID,1)-mpPos(mp,1))*(vtxCoords(vID,2)-mpPos(mp,2))/(radius*radius));
          Kokkos::atomic_add(&m44_d(i,elm), w_vtx*mScale*(vtxCoords(vID,2)-mpPos(mp,2))*(vtxCoords(vID,2)-mpPos(mp,2))/(radius*radius));
      }
    }
  };
  p_MPs->parallel_for(sub_assemble, "sub_assembly");
  pumipic::RecordTime("VtxSubAssemblyComputeCoeff", timer.seconds());

  Kokkos::Timer timer2;
  kkDbl2dViewHostU m11_h(m11, vtxPerElm, nCells);
  kkDbl2dViewHostU m12_h(m12, vtxPerElm, nCells);
  kkDbl2dViewHostU m13_h(m13, vtxPerElm, nCells);
  kkDbl2dViewHostU m14_h(m14, vtxPerElm, nCells);
  kkDbl2dViewHostU m22_h(m22, vtxPerElm, nCells);
  kkDbl2dViewHostU m23_h(m23, vtxPerElm, nCells);
  kkDbl2dViewHostU m24_h(m24, vtxPerElm, nCells);
  kkDbl2dViewHostU m33_h(m33, vtxPerElm, nCells);
  kkDbl2dViewHostU m34_h(m34, vtxPerElm, nCells);
  kkDbl2dViewHostU m44_h(m44, vtxPerElm, nCells);
  
  Kokkos::deep_copy(m11_h, m11_d); 
  Kokkos::deep_copy(m12_h, m12_d); 
  Kokkos::deep_copy(m13_h, m13_d); 
  Kokkos::deep_copy(m14_h, m14_d); 
  Kokkos::deep_copy(m22_h, m22_d); 
  Kokkos::deep_copy(m23_h, m23_d); 
  Kokkos::deep_copy(m24_h, m24_d); 
  Kokkos::deep_copy(m33_h, m33_d); 
  Kokkos::deep_copy(m34_h, m34_d); 
  Kokkos::deep_copy(m44_h, m44_d); 
  pumipic::RecordTime("VtxSubAssemblyGetCoeff", timer2.seconds());
  
}

void MPMesh::solveMatrixAndRegularize(int nVertices, double* m11, double* m12, double* m13, double* m14, 
                                       double* m22, double* m23, double* m24, 
                                       double* m33, double* m34,
                                       double* m44){

  Kokkos::Timer timer;
  kkViewHostU<const double*> m11_h(m11, nVertices);
  kkViewHostU<const double*> m12_h(m12, nVertices);
  kkViewHostU<const double*> m13_h(m13, nVertices);
  kkViewHostU<const double*> m14_h(m14, nVertices);
  kkViewHostU<const double*> m22_h(m22, nVertices);
  kkViewHostU<const double*> m23_h(m23, nVertices);
  kkViewHostU<const double*> m24_h(m24, nVertices);
  kkViewHostU<const double*> m33_h(m33, nVertices);
  kkViewHostU<const double*> m34_h(m34, nVertices);
  kkViewHostU<const double*> m44_h(m44, nVertices);
  
  Kokkos::View<double*> m11_d("m11", nVertices);
  Kokkos::View<double*> m12_d("m12", nVertices);
  Kokkos::View<double*> m13_d("m13", nVertices);
  Kokkos::View<double*> m14_d("m14", nVertices); 
  Kokkos::View<double*> m22_d("m22", nVertices);
  Kokkos::View<double*> m23_d("m23", nVertices);
  Kokkos::View<double*> m24_d("m24", nVertices);
  Kokkos::View<double*> m33_d("m33", nVertices);
  Kokkos::View<double*> m34_d("m34", nVertices);
  Kokkos::View<double*> m44_d("m44", nVertices);
  
  Kokkos::deep_copy(m11_d, m11_h);
  Kokkos::deep_copy(m12_d, m12_h);
  Kokkos::deep_copy(m13_d, m13_h);
  Kokkos::deep_copy(m14_d, m14_h);
  Kokkos::deep_copy(m22_d, m22_h);
  Kokkos::deep_copy(m23_d, m23_h);
  Kokkos::deep_copy(m24_d, m24_h);
  Kokkos::deep_copy(m33_d, m33_h);
  Kokkos::deep_copy(m34_d, m34_h);
  Kokkos::deep_copy(m44_d, m44_h);
  pumipic::RecordTime("polyMPOsolveMatrixCoeffSet", timer.seconds());

  Kokkos::Timer timer2;
  auto dual_triangle_area=p_mesh->getMeshField<MeshF_DualTriangleArea>();
  Kokkos::View<double*[vec4d_nEntries]> VtxCoeffs("VtxCoeffs", nVertices);
  double radius=p_mesh->getSphereRadius();
  Kokkos::parallel_for("fill", nVertices, KOKKOS_LAMBDA(const int vtx){
    Vec4d v0 = {m11_d(vtx), m12_d(vtx), m13_d(vtx), m14_d(vtx)};
    Vec4d v1 = {m12_d(vtx), m22_d(vtx), m23_d(vtx), m24_d(vtx)};
    Vec4d v2 = {m13_d(vtx), m23_d(vtx), m33_d(vtx), m34_d(vtx)};
    Vec4d v3 = {m14_d(vtx), m24_d(vtx), m34_d(vtx), m44_d(vtx)}; 
    //Matrix4d A = {v0,v1,v2,v3};
    Matrix4d A_regularized = {v0, v1, v2, v3};
    double regParam = sqrt(EPSILON)*(m11_d(vtx) + m22_d(vtx) + m33_d(vtx) + m44_d(vtx));
    A_regularized.addToDiag(regParam);
    
    double coeff[vec4d_nEntries]={0.0, 0.0, 0.0, 0.0};
    CholeskySolve4d_UnitRHS(A_regularized, coeff);
    
    double mScale=sqrt(dual_triangle_area(vtx,0))/radius;
    coeff[0]=coeff[0]*mScale*mScale;
    coeff[1]=coeff[1]*mScale;
    coeff[2]=coeff[2]*mScale;
    coeff[3]=coeff[3]*mScale;

    for (int i=0; i<vec4d_nEntries; i++) 
      VtxCoeffs(vtx,i)=coeff[i];
  });
  this->precomputedVtxCoeffs = VtxCoeffs;
  pumipic::RecordTime("polyMPOsolveMatrixCoeffCompute", timer2.seconds());

}

template <MeshFieldIndex meshFieldIndex>
void MPMesh::subAssemblyVtx1(int vtxPerElm, int nCells, int comp, double* array) {
  Kokkos::Timer timer; 
  
  auto VtxCoeffs=this->precomputedVtxCoeffs; 
  
  // MPI Information
  MPI_Comm comm = p_MPs->getMPIComm(); 
  int comm_rank;
  MPI_Comm_rank(comm, &comm_rank);

  //Mesh Information
  auto elm2VtxConn = p_mesh->getElm2VtxConn();  
  int numVtx = p_mesh->getNumVertices();
  auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
  auto elm2Process = p_mesh->getElm2Process();
 
  // Material Points Information
  constexpr MaterialPointSlice mpfIndex = meshFieldIndexToMPSlice<meshFieldIndex>;
  auto mpData = p_MPs->getData<mpfIndex>();
  auto weight = p_MPs->getData<MPF_Basis_Vals>();
  auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
 
  double radius=p_mesh->getSphereRadius();

  Kokkos::View<double**> array_d("reconstructedField", vtxPerElm, nCells);
  auto sub_assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask && (elm2Process(elm)==comm_rank)) { 
      int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
      for(int i=0; i<nVtxE; i++){
        int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
        double w_vtx=weight(mp,i);
        double CoordDiffs[vec4d_nEntries] = {1, (vtxCoords(vID,0) - mpPositions(mp,0))/radius, 
                                                (vtxCoords(vID,1) - mpPositions(mp,1))/radius, 
                                                (vtxCoords(vID,2) - mpPositions(mp,2))/radius};

        auto factor = w_vtx*(VtxCoeffs(vID,0) + VtxCoeffs(vID,1)*CoordDiffs[1] + 
                                                VtxCoeffs(vID,2)*CoordDiffs[2] + 
                                                VtxCoeffs(vID,3)*CoordDiffs[3]);
  
        auto val = factor*mpData(mp, comp);
        Kokkos::atomic_add(&array_d(i, elm), val);
      }
    }
  };
  p_MPs->parallel_for(sub_assemble, "sub_assembly"); 
  pumipic::RecordTime("polyMPOsubAssemblyFieldCompute", timer.seconds());
  
  Kokkos::Timer timer2;
  kkDbl2dViewHostU arrayHost(array, vtxPerElm, nCells);
  Kokkos::deep_copy(arrayHost, array_d); 
  pumipic::RecordTime("PolyMPOsubAssemblyFieldGet", timer2.seconds());
}

// An improvement on the above method by doing the full assembly on GPUs
void MPMesh::assembleField(int vtxPerElm, int nCells, int nVerticesSolve, int nVertices, double* array_sub, double* array_full){
  
  //Mesh Information
  auto elm2VtxConn = p_mesh->getElm2VtxConn();  
  int numVtx = p_mesh->getNumVertices();
  auto elm2Process = p_mesh->getElm2Process();

  //Copy the subAssembled Field to GPU
  Kokkos::Timer timer;
  kkViewHostU<const double**> array_sub_h(array_sub, vtxPerElm, nCells);
  Kokkos::View<double**> array_sub_d("array_sub", vtxPerElm, nCells);
  Kokkos::deep_copy(array_sub_d, array_sub_h);
  pumipic::RecordTime("polyMPOsetsubAssemblyField", timer.seconds());

  Kokkos::Timer timer1;
  Kokkos::View<double*> array_full_d("reconstructedField", nVertices);
  Kokkos::parallel_for("assemble", nCells, KOKKOS_LAMBDA(const int elm){
    int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
    for(int i=0; i<nVtxE; i++){
      int vID = elm2VtxConn(elm,i+1)-1;
      if(vID < nVerticesSolve){
        auto val = array_sub_d(i, elm);
        Kokkos::atomic_add(&array_full_d(vID), val);
      }
    } 
  });
  pumipic::RecordTime("polyMPOfullAssemble", timer1.seconds());

  //Copy the assembled field to CPU
  Kokkos::Timer timer2;
  kkDblViewHostU arrayHost(array_full, nVertices);
  Kokkos::deep_copy(arrayHost, array_full_d);
  pumipic::RecordTime("polyMPOgetAssemblyField", timer2.seconds());

  MPMesh::startCommunication(nVerticesSolve);
}

//Start Communication routine
void MPMesh::startCommunication(int nVerticesSolve){

  MPI_Comm comm = p_MPs->getMPIComm(); 
  int comm_rank, nProcs;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &nProcs); 

  int nCells = p_mesh->getNumElements();
  int numVertices = p_mesh->getNumVertices();

  //Owning processes in GPU and copy to CPU
  auto elmOwners = p_mesh->getElm2Process();
  auto elmOwners_host = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace::memory_space(),
                        elmOwners);

  auto vtxGlobal= p_mesh->getVtxGlobal();
  auto vtxGlobal_host = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace::memory_space(),
                        vtxGlobal);

  //Elment to Vertex Connection in GPU and communicate to CPU
  auto elm2VtxConn = p_mesh->getElm2VtxConn();  
  auto elm2VtxConn_host = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace::memory_space(),
                          elm2VtxConn);
  //Debugging
  for (auto k=0; k<elmOwners_host.size(); k++)
    if(k==0 || k==elmOwners_host.size()-1)
      printf("Rank %d Owning Proc %d \n", comm_rank, elmOwners_host(k));
  
  //Find adjacent processes and no of adjacent cells for each process
  std::map<int, int> adjProcsForCells;
  std::vector<int> send_each_proc(nProcs);
  for (auto iCell=0; iCell<nCells; iCell++){
    if(elmOwners_host(iCell) != comm_rank){
      auto ownerProc = elmOwners_host[iCell];
      adjProcsForCells[ownerProc] = adjProcsForCells[ownerProc]+1;
      send_each_proc[ownerProc] +=1;
    }
  }
  //Debugging
  printf("Size adjProcs %d \n", adjProcsForCells.size());
  for (const auto& [key, value] : adjProcsForCells){
    printf ("Process %d sends to %d size %d\n", comm_rank, key, value);
  }
 
  //Find received particles in each process
  std::vector<int> total_recv_each_proc(nProcs);
  MPI_Alltoall(send_each_proc.data(), 1, MPI_INT, total_recv_each_proc.data(), 1, MPI_INT, comm);
  //Debiugging for receiving
  for (int i=0; i<nProcs; i++)
    printf("Rank %d receiving from %d size %d \n", comm_rank, i, total_recv_each_proc[i]);

  //Create maps of data to Send
  std::map<int, std::vector<int>> cellDataToSend;
  std::map<int, int> counter;
  for (const auto& [key, value] : adjProcsForCells){
    cellDataToSend[key].resize(4 * value * maxVtxsPerElm);
    counter[key] = 0;
  }
  
  for (auto iCell=0; iCell<nCells; iCell++){
    if(elmOwners_host(iCell) != comm_rank){
      int ownerProc = elmOwners_host[iCell];
      auto idx_start = counter[ownerProc]*4*maxVtxsPerElm;
      int nVtxE = elm2VtxConn_host(iCell,0);
      for (int v=0; v<nVtxE; v++){
        int vID = elm2VtxConn_host(iCell, v+1)-1;
        int idx = idx_start + v*4;
        if (vID < nVerticesSolve)
          cellDataToSend[ownerProc][idx+0] = 0;            //TO DO better way
        else
          cellDataToSend[ownerProc][idx+0] = 1;
        cellDataToSend[ownerProc][idx+1] = comm_rank;      //sending Proc TODO not needed
        cellDataToSend[ownerProc][idx+2] = vID;            //localID
        cellDataToSend[ownerProc][idx+3] = vtxGlobal_host(vID); //globalID
      }
      counter[ownerProc] = counter[ownerProc] + 1;
      //assert(counter[ownerProc] == adjProcsForCells.find(ownerProc));
    }
  }
  
  std::vector<MPI_Request> s_requests;
  s_requests.resize(cellDataToSend.size());
  int count_s_request=0;
  for (auto & [proc, vec] : cellDataToSend){
    MPI_Isend(vec.data(), vec.size(), MPI_INT, proc, MPI_ANY_TAG, comm, &s_requests[count_s_request]);
    count_s_request=count_s_request+1;
  }

  std::vector<std::vector<int>> cellDataToReceive;
  cellDataToReceive.resize(nProcs); //
  for (int iProc=0; iProc< nProcs; iProc++)
    cellDataToReceive[iProc].resize(total_recv_each_proc[iProc]);   
  
  std::vector<MPI_Request> r_requests;
  r_requests.resize(nProcs);
  for (int iProc=0; iProc< nProcs; iProc++)
    MPI_Irecv(cellDataToReceive[iProc].data(), total_recv_each_proc[iProc], MPI_INT, iProc, 
              MPI_ANY_TAG, comm, &r_requests[iProc]);
  
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
