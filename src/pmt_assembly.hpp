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
    
    DoubleView vField("vField2",numVtxs);
    auto MPs = mpMesh.MPs;
    auto mpPositions = MPs->getData<0>();
    auto assemble = PS_LAMBDA(const int& e, const int& pid, const int& mask) {
      if(mask > 0) {
        int nVtxE = elm2VtxConn(e,0);
        for(int i=0; i<nVtxE; i++){
            int vID = elm2VtxConn(e,i+1)-1;
            auto vertexLoc = vtxCoords(vID);
            double distance = mpPositions(pid,0);
            Kokkos::atomic_add(&vField(vID),distance);
        }
      }
    };
    MPs->parallel_for(assemble);

   
    return vField;
}

} //end namespace polyMpmTest
#endif
