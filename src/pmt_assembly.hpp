#ifndef POLYMPMTEST_ASSEMBLY_H
#define POLYMPMTEST_ASSEMBLY_H

#include "pmt_MPMesh.hpp"
#include "pmt_wachspressBasis.hpp"

namespace polyMpmTest{

DoubleView assemblyV0(MPMesh& mpMesh){
    auto mesh = *mpMesh.mesh;
    int numVtxs = mesh.getNumVertices();
    auto elm2VtxConn = mesh.getElm2VtxConn();
    
    DoubleView vField("vField2",numVtxs);
    auto MPs = mpMesh.MPs;
    auto mpPositions = MPs->getData<MPF_Cur_Pos_XYZ>(); //get the array of MP coordinates/positions
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
    MPs->parallel_for(assemble, "assembly");
    return vField;
}

template <MaterialPointSlice mpfIndex, MeshFieldIndex mfIndex>
void assembly(MPMesh& mpMesh, bool basisWeightFlag, bool massWeightFlag){
    if(basisWeightFlag || massWeightFlag) {
      std::cerr << "WARNING: basis and mass weight flags ignored\n";
    }
    auto mesh = *mpMesh.mesh;
    int numVtxs = mesh.getNumVertices();
    auto elm2VtxConn = mesh.getElm2VtxConn();
   
    auto MPs = mpMesh.MPs;
    auto mpData = MPs->getData<mpfIndex>();
    auto massWeight = MPs->getData<MPF_Mass>();
    //TODO:massWeight is not used in the loop
    //if(!massWeightFlag){
        //massWeight =  
    //}
    auto basis = MPs->getData<MPF_Basis_Vals>();
    //if(!basisWeightFlag){
    //}
    const int numEntries = mpSlice2MeshFieldIndex.at(mpfIndex).first;
    const MeshFieldIndex meshFieldIndex = mpSlice2MeshFieldIndex.at(mpfIndex).second;
    PMT_ALWAYS_ASSERT(meshFieldIndex == mfIndex);
    auto meshField = mesh.getMeshField<mfIndex>(); 
    //auto meshField = mesh.getMeshField<Mesh_Field_Cur_Pos_XYZ>(); 
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
    MPs->parallel_for(assemble, "assembly");
}

// (HDT) weighted assembly of scalar field
template <MaterialPointSlice index>
DoubleView wtScaAssembly(MPMesh& mpMesh){
    auto mesh = *mpMesh.mesh;
    auto vtxCoords = mesh.getVtxCoords();
    int numVtxs = mesh.getNumVertices(); // total number of vertices of the mesh
    auto elm2VtxConn = mesh.getElm2VtxConn();
    auto MPs = mpMesh.MPs;
    auto mpPositions = MPs->getData<MPF_Cur_Pos_XYZ>();
    
    DoubleView vField("wtScaField", numVtxs); // Kokkos array of double type, size = numVtxs

    auto mpData = MPs->getData<index>();
    
    auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
      if (mask) {
        /* get the coordinates of all the vertices of elm */
        int nElmVtxs = elm2VtxConn(elm,0);      // number of vertices bounding the element
        Vector2 eVtxCoords[maxVtxsPerElm + 1];
        for (int i = 1; i <= nElmVtxs; i++) {
          eVtxCoords[i-1] = vtxCoords(elm2VtxConn(elm,i)-1);    // elm2VtxConn(elm,i) is the vertex ID (1-based index) of vertex #i of elm
        }
        eVtxCoords[nElmVtxs] = vtxCoords(elm2VtxConn(elm,1)-1); // last component of eVtxCoords stores the firs vertex (to avoid if-condition in the Wachspress computation)
        
        /* compute the values of basis functions at mp position */
        double basisByArea[maxElmsPerVtx];
        Vector2 mpCoord(mpPositions(mp,0), mpPositions(mp,1));
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
    MPs->parallel_for(assemble, "weightedScalarAssembly");
    return vField;
} // wtScaAssembly


// (HDT) weighted assembly of vector2 field (not weighted by mass/volume)
template <MaterialPointSlice index>
Vector2View wtVec2Assembly(MPMesh& mpMesh){
    auto mesh = *mpMesh.mesh;
    auto vtxCoords = mesh.getVtxCoords();
    int numVtxs = mesh.getNumVertices(); // total number of vertices of the mesh
    auto elm2VtxConn = mesh.getElm2VtxConn();
    auto MPs = mpMesh.MPs;
    auto mpPositions = MPs->getData<MPF_Cur_Pos_XYZ>();
    
    Vector2View vField("wtVec2Field", numVtxs); // Kokkos array of Vector2 type, size = numVtxs

    auto mpData = MPs->getData<index>();
    
    auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
      if (mask) {
        /* collect the coordinates of all the vertices of elm */
        int nElmVtxs = elm2VtxConn(elm,0);      // number of vertices bounding the element
        Vector2 eVtxCoords[maxVtxsPerElm + 1];
        for (int i = 1; i <= nElmVtxs; i++) {
          eVtxCoords[i-1] = vtxCoords(elm2VtxConn(elm,i)-1);    // elm2VtxConn(elm,i) is the vertex ID (1-based index) of vertex #i of elm
        }
        eVtxCoords[nElmVtxs] = vtxCoords(elm2VtxConn(elm,1)-1); // last component of eVtxCoords stores the firs vertex (to avoid if-condition in the Wachspress computation)
        
        /* compute the values of basis functions at mp position */
        double basisByArea[maxElmsPerVtx];
        Vector2 mpCoord(mpPositions(mp,0), mpPositions(mp,1));
        getBasisByAreaGblForm(mpCoord, nElmVtxs, eVtxCoords, basisByArea);

        /* get the mp's volume */
        double mpVolume = 1.0; // TODO: change to mp's volume here

        /* get the mp's property to be assembled */
        Vector2 assValue;
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
    MPs->parallel_for(assemble, "weightedVector2Assembly");
    return vField;
} // wtVec2Assembly

} //end namespace polyMpmTest
#endif
