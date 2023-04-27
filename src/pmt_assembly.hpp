#ifndef POLYMPMTEST_ASSEMBLY_H
#define POLYMPMTEST_ASSEMBLY_H

#include "pmt_MPMesh.hpp"

namespace polyMpmTest{

DoubleView assembly(MPMesh& mpMesh){
    auto mesh = mpMesh.getMesh();
    int numVtxs = mesh.getNumVertices();
    auto elm2VtxConn = mesh.getElm2VtxConn();
    
    DoubleView vField("vField2",numVtxs);
    auto MPs = mpMesh.MPs;
    auto mpPositions = MPs->getData<0>();
    auto assemble = PS_LAMBDA(const int& elm, const int& pid, const int& mask) {
      if(mask > 0) {
        int nVtxE = elm2VtxConn(elm,0);
        for(int i=0; i<nVtxE; i++){
          int vID = elm2VtxConn(elm,i+1)-1;
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
