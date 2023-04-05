#include <Kokkos_Core.hpp>
#include "pmt_utils.hpp"
#include "pmt_MPM.hpp"

namespace polyMpmTest{

const bool printVTP = true;

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
   
    //numMPs = 10;//XXX
    Vector2View history("positionHistory",numMPs);
    Vector2View::HostMirror h_history = Kokkos::create_mirror_view(history);
    Vector2View resultLeft("positionResult",numMPs);
    Vector2View::HostMirror h_resultLeft = Kokkos::create_mirror_view(resultLeft);
    Vector2View resultRight("positionResult",numMPs);
    Vector2View::HostMirror h_resultRight = Kokkos::create_mirror_view(resultRight);
    Vector2View::HostMirror h_MPsPosition = Kokkos::create_mirror_view(MPsPosition);
    Kokkos::parallel_for("test",numMPs,KOKKOS_LAMBDA(const int iMP){
        Vector2 MP = MPsPosition(iMP);
        if(isActive(iMP)){
            int iElm = MPs2Elm(iMP);
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
                        //printf("%d: MP=(%.3f,%.3f) MPnew=(%.3f,%.3f) MATLAB: [%f %f], [%f %f]\n",iMP, MP[0], MP[1], MPnew[0], MPnew[1], MP[0], MPnew[0], MP[1], MPnew[1]);
                        //go to the next elm
                        //int iElmOld = iElm;
                        iElm = elm2ElmConn(iElm,i+1);
                        //if(MP[0]-464621<1 && MP[0]-464621>0){
                        //    printf("%d: %f*%f= %f && eiCross = %f\n",i,pdx[i],pdx[ip1], pdx[i]*pdx[ip1] , e[i].cross(MPnew-vtxCoords(v[i])));
                        //    printf("%d: from %d to %d, MP= (%f,%f), dx= (%f,%f)\n",iMP ,iElmOld, iElm , MP[0], MP[1], dx(iMP)[0], dx(iMP)[1]);
                        //}
                        //printf("%d: from %d to %d, MP= (%f,%f), dx= (%f,%f)\n",iMP ,iElmOld, iElm , MP[0], MP[1], dx(iMP)[0], dx(iMP)[1]);
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
                if(printVTP){ 
                    Vector2 MParrow = MP + dx(iMP)*0.7;
                    Vector2 shift = Vector2(-dx(iMP)[1],dx(iMP)[0])*0.1;
                    Vector2 MPLeft = MParrow + shift;
                    Vector2 MPRight = MParrow - shift;
                    Kokkos::atomic_store(&history(iMP),MP);
                    Kokkos::atomic_store(&resultLeft(iMP),MPLeft);
                    Kokkos::atomic_store(&resultRight(iMP),MPRight);
                }
                MPs2Elm(iMP) = iElm;
                MPsPosition(iMP) = MPnew;
               break;
            }
        }
        else{
            Kokkos::atomic_store(&history(iMP),MP);
            Kokkos::atomic_store(&resultLeft(iMP),MP);
            Kokkos::atomic_store(&resultRight(iMP),MP);
        }
    }); 
    if(printVTP){ 
        printf("<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"%d\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",numMPs*4,numMPs*2); 
        Kokkos::fence();
        Kokkos::deep_copy(h_history, history);
        Kokkos::deep_copy(h_resultLeft, resultLeft);
        Kokkos::deep_copy(h_resultRight, resultRight);
        Kokkos::deep_copy(h_MPsPosition, MPsPosition);
        Kokkos::fence();
        for(int i=0; i<numMPs; i++){
            //XXX: MPsPosition is the updated new position, h_history is the old position
            printf("          %f %f 0.0\n          %f %f 0.0\n          %f %f 0.0\n          %f %f 0.0\n",h_history(i)[0], h_history(i)[1], h_MPsPosition(i)[0], h_MPsPosition(i)[1], h_resultLeft(i)[0], h_resultLeft(i)[1], h_resultRight(i)[0], h_resultRight(i)[1]);
        }
        printf("        </DataArray>\n      </Points>\n      <Lines>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"); 
        for(int i=0; i<numMPs*4; i+=4){
            // 01 213
            printf("          %d %d\n          %d %d %d\n", i, i+1, i+2, i+1, i+3);
        }
        printf("        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        for(int i=0; i<numMPs*5; i+=5){
            printf("          %d\n          %d\n",i+2,i+5);
        }
        printf("        </DataArray>\n      </Lines>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");

   }

}

} 
