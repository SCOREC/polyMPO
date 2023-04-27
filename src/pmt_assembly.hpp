#ifndef POLYMPMTEST_ASSEMBLY_H
#define POLYMPMTEST_ASSEMBLY_H

#include "pmt_MPMesh.hpp"

namespace polyMpmTest{

DoubleView assembly(MPMesh& mpMesh){
    auto mesh = mpMesh.getMesh();
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    int numMPs = mpMesh.MPs->getCount();
     
    auto vtxCoords = mesh.getVtxCoords(); 
    auto elm2VtxConn = mesh.getElm2VtxConn();
    auto xp = mpMesh.MPs->getPositions();
    
    auto MPs2Elm = mpMesh.getMPs2Elm();
    DoubleView vField("vField2",numVtxs);
    Kokkos::parallel_for("vertex_assem", numMPs, KOKKOS_LAMBDA(const int iMP){
        int ielm = MPs2Elm(iMP); 
        int nVtxE = elm2VtxConn(ielm,0);
        for(int i=0; i<nVtxE; i++){
            int vID = elm2VtxConn(ielm,i+1)-1;
            auto vertexLoc = vtxCoords(vID);
            double distance = xp(iMP)[0];
            Kokkos::atomic_add(&vField(vID),distance);
        }
    });
   
    return vField;
}

} //end namespace polyMpmTest
#endif
