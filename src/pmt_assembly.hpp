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
    auto mpPositions = MPs->getData<MP_CUR_POS_XYZ>(); //get the array of MP coordinates/positions
    auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    //for elm in elementsInMesh { //pseudo code - the 'parallel_for' handles this
    //  for mp in materialPointsInElm { //pseudo code (cont.)
          if(mask) { //if material point is 'active'/'enabled'
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

template <MaterialPointSlice index>
DoubleView assemblyNew(MPMesh& mpMesh){
    auto mesh = mpMesh.getMesh();
    int numVtxs = mesh.getNumVertices();
    auto elm2VtxConn = mesh.getElm2VtxConn();
    
    DoubleView vField("vField",numVtxs);
    auto MPs = mpMesh.MPs;
    auto mpData = MPs->getData<index>();
    int loopNum = 0;
    switch(index){
        case MP_CUR_POS_LAT_LON:
        case MP_TGT_POS_LAT_LON:
        case MP_VEL:
            loopNum = 2;
            break;
        case MP_CUR_POS_XYZ:
        case MP_TGT_POS_XYZ:
        case MP_STRESS_DIV:
        case MP_SHEAR_TRACTION:
            loopNum = 3;
            break;
        case MP_STRAIN_RATE:
        case MP_STRESS:
            loopNum = 6;
            break;
        case MP_BASIS_VALS:
            loopNum = 8;
            break;
        case MP_CONSTV_MDL_PARAM:
            loopNum = 12;
            break;
        //these case won't be assemble
        case MP_STATUS:
        case MP_CUR_ELM_ID:
        case MP_TGT_ELM_ID:
        case MP_FLAG_BASIS_VALS:
        case MP_BASIS_GRAD_VALS:
            PMT_ALWAYS_ASSERT(false);
    }   
    auto assemble = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
        if(mask) { //if material point is 'active'/'enabled'
            int nVtxE = elm2VtxConn(elm,0); //number of vertices bounding the element
            for(int i=0; i<nVtxE; i++){
                int vID = elm2VtxConn(elm,i+1)-1; //vID = vertex id
                double distance = 0;
                for(int i=0;i<loopNum;i++)
                    distance += mpData(mp,0);
                Kokkos::atomic_add(&vField(vID),distance);
            }
          }
    };
    MPs->parallel_for(assemble, "assembly");
    return vField;
}

} //end namespace polyMpmTest
#endif
