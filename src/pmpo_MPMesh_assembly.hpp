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

template <MaterialPointSlice mpfIndex>
void MPMesh::assembly(bool basisWeightFlag, bool massWeightFlag){
    if(basisWeightFlag || massWeightFlag) {
      std::cerr << "WARNING: basis and mass weight flags ignored\n";
    }
    auto elm2VtxConn = p_mesh->getElm2VtxConn();
   
    auto mpData = p_MPs->getData<mpfIndex>();
    auto massWeight = p_MPs->getData<MPF_Mass>();
    //TODO:massWeight is not used in the loop
    //if(!massWeightFlag){
        //massWeight =  
    //}
    auto basis = p_MPs->getData<MPF_Basis_Vals>();
    //if(!basisWeightFlag){
    //}
    const int numEntries = mpSliceToMeshFieldSize<mpfIndex>;
    constexpr MeshFieldIndex meshFieldIndex = mpSliceToMeshFieldIndex<mpfIndex>;
    auto meshField = p_mesh->getMeshField<meshFieldIndex>(); 
    //auto meshField = p_mesh->getMeshField<Mesh_Field_Cur_Pos_XYZ>(); 
    auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if(mask) { //if material point is 'active'/'enabled'
            int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
            double mpMass = massWeight(mp,0);
            for(int i=0; i<nVtxE; i++){
                int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
                double fieldComponentVal;
                for(int j=0;j<numEntries;j++){
                    fieldComponentVal = mpData(mp,j)*basis(mp,i)*mpMass;
                    Kokkos::atomic_add(&meshField(vID,j),fieldComponentVal);
                }
            }
        }
    };
    p_MPs->parallel_for(assemble, "assembly");
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

template<MaterialPointSlice mpSliceIndex>
void MPMesh::setReconstructSlice() {
  auto function = [this](){ assembly<mpSliceIndex>(false, false); }; 
  const auto [iter, success] = reconstructSlice.insert({mpSliceIndex, function});
  if (!success){
    std::cerr << "Error: Slice is already being reconstructed\n";
    exit(1);
  }
}

} //end namespace polyMPO
#endif
