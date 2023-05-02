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
    auto mpPositions = MPs->getData<0>(); //get the array of MP coordinates/positions
    auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    //for elm in elementsInMesh { //pseudo code - the 'parallel_for' handles this
    //  for mp in materialPointsInElm { //pseudo code (cont.)
          if(mask > 0) { //if material point is 'active'/'enabled'
            int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
            for(int i=0; i<nVtxE; i++){
              int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
              double distance = mpPositions(mp,0);
              Kokkos::atomic_add(&vField(vID),distance);
            }
          }
    //  }
    //}
    };
    MPs->parallel_for(assemble, "assembly");
    return vField;
}

} //end namespace polyMpmTest
#endif
