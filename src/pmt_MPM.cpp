#include <Kokkos_Core.hpp>
#include "pmt_utils.hpp"
#include "pmt_MPM.hpp"

namespace polyMpmTest{

void MPM::T2LTracking(Vector2View dx){
    int numVtxs = mesh_.getNumVertices();
    int numElms = mesh_.getNumElements();
    
    const auto vtxCoords = mesh_.getVtxCoords(); 
    auto elm2VtxConn = mesh_.getElm2VtxConn();
    auto vtx2ElmConn = mesh_.getVtx2ElmConn();
    auto elm2ElmConn = mesh_.getElm2ElmConn();

    auto numMPs = materialPoints_.getCount();
    auto MPsPosition = materialPoints_.getPositions();
    auto isActive = materialPoints_.isActive();
    
    auto MPs2Elm = materialPoints2Elm_;
   
    //numMPs = 1;//XXX 
    Kokkos::parallel_for("test",numMPs,KOKKOS_LAMBDA(const int iMP){
        if(isActive(iMP)){
            int iElm = MPs2Elm(iMP);
            Vector2 MP = MPsPosition(iMP);
            Vector2 MPnew = MP + dx(iMP);    
            //printf(" MP=(%.3f,%.3f) MPnew=(%.3f,%.3f) MATLAB: [%f %f], [%f %f]\n", MP[0], MP[1], MPnew[0], MPnew[1], MP[0], MPnew[0], MP[1], MPnew[1]);

            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                bool goToNeighbour = false;
                //seperate the elm2Vtx
                int v[maxVtxsPerElm];
                for(int i=0; i< numVtx; i++)
                    v[i] = elm2VtxConn(iElm,i+1)-1;
                //get edges and perpendiculardx
                Vector2 e[maxVtxsPerElm];
                double pdx[maxVtxsPerElm];                    
                for(int i=0; i< numVtx; i++){
                    Vector2 v_i = vtxCoords(v[i]);
                    Vector2 v_ip1 = vtxCoords(v[(i+1)%numVtx]);
                    e[i] = v_ip1 - v_i;
                    pdx[i] = (v_i - MP).cross(dx(iMP));
                }
                
                for(int i=0; i<numVtx; i++){
                    int ip1 = (i+1)%numVtx;
                    //pdx*pdx<0 and
                    if(pdx[i]*pdx[ip1] <0 && e[i].cross(MPnew-vtxCoords(v[i]))<0){
                        printf("%d: MP=(%.3f,%.3f) MPnew=(%.3f,%.3f) MATLAB: [%f %f], [%f %f]\n",iMP, MP[0], MP[1], MPnew[0], MPnew[1], MP[0], MPnew[0], MP[1], MPnew[1]);
                        //go to the next elm
                        int iElmOld = iElm;
                        iElm = elm2ElmConn(iElm,i+1);
                        printf("%d: from %d to %d, MP= (%f,%f), dx= (%f,%f)\n",iMP ,iElmOld, iElm , MP[0], MP[1], dx(iMP)[0], dx(iMP)[1]);
                        goToNeighbour = true;
                        if(iElm <0){
                            isActive(iMP) = false;
                            goToNeighbour = false;
                        }
                    }
                }
                //if goes to the other 
                if(goToNeighbour)
                    continue; 
                //otherwise we do the update and end the loop
                MPs2Elm(iMP) = iElm;
                MPsPosition(iMP) = MPnew;
                break;
            }
        }
        //printf("position update: (%f,%f)\n",MPsPosition(0)[0], MPsPosition(0)[1]);
    });   
} 

} 
