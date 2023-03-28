#include <Kokkos_Core.hpp>
#include "pmt_utils.hpp"
#include "pmt_MPM.hpp"

namespace polyMpmTest{

void MPM::T2LTracking(Vector2View dx){
    int numVtxs = mesh_.getNumVertices();
    int numElms = mesh_.getNumElements();
    
    auto vtxCoords = mesh_.getVtxCoords(); 
    auto elm2VtxConn = mesh_.getElm2VtxConn();
    auto vtx2ElmConn = mesh_.getVtx2ElmConn();

    auto numMPs = materialPoints_.getCount();
    auto MPsPosition = materialPoints_.getPositions();
    auto isActive = materialPoints_.isActive();
    
    Kokkos::parallel_for("test",numMPs,KOKKOS_LAMBDA(const int i){
        if (isActive(iMP)){
            int iElm = MPs2Elm(iMP);
            Vector2 v[maxVtxsPerElm+1] = {vtxCoords(elm2VtxConn(iElm,1))};
            initArray(v,maxVtxsPerElm+1,vtxCoords(elm2VtxConn(iElm,1)));
            int numVtx = elm2VtxConn(iElm,0);
            for(int i = 1; i<=numVtx; i++){
                v[i-1] = vtxCoords(elm2VtxConn(iElm,i)-1);
            }
            v[numVtx] = vtxCoords(elm2VtxConn(iElm,1)-1);
                
            //get edges
        }
    });   
} 

} 
