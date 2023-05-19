#include <Kokkos_Core.hpp>
#include "pmt_utils.hpp"
#include "pmt_MPM.hpp"

namespace polyMpmTest{

void MPM::CVTTrackingEdgeCenterBased(Vector2View dx, const int printVTP){
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

    //numMPs = 20; //XXX
    IntView count("countCrossMPs",numMPs);
    Kokkos::View<Vector2*[maxVtxsPerElm]> edgeCenters("EdgeCenters",numElms);
    Kokkos::parallel_for("EdgeCenterCalc",numElms,KOKKOS_LAMBDA(const int iElm){
        int numVtx = elm2VtxConn(iElm,0);
        int v[maxVtxsPerElm];
        for(int i=0; i< numVtx; i++)
            v[i] = elm2VtxConn(iElm,i+1)-1;
        for(int i=0; i< numVtx; i++){
            Vector2 v_i = vtxCoords(v[i]);
            Vector2 v_ip1 = vtxCoords(v[(i+1)%numVtx]);
            edgeCenters(iElm,i) = (v_ip1 + v_i)*0.5;
        }
   });
    Kokkos::fence();
    Kokkos::parallel_for("CVTEdgeCalc",numMPs,KOKKOS_LAMBDA(const int iMP){
        Vector2 MP = MPsPosition(iMP);
        if(isActive(iMP)){
            int iElm = MPs2Elm(iMP);
            Vector2 MPnew = MP + dx(iMP);
            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                //calc dist square from each edge center to MPnew
                //calc dot products to check inside or not
                int edgeIndex = -1;
                double minDistSq = DBL_MAX;
                for(int i=0; i< numVtx; i++){
                    Vector2 edgeCenter = edgeCenters(iElm,i);
                    Vector2 delta = MPnew - edgeCenter;
                    double currentDistSq = delta[0]*delta[0] + delta[1]*delta[1];
                    double dotProduct = dx(iMP).dot(delta);
                        if(iMP == 29){
                        Vector2 c = calcElmCenter(iElm,elm2VtxConn,vtxCoords);
                        //printf("iMP iElm:%d %d %f %f\n",iMP,iElm,c[0],c[1]);
                        //printf("iMP iElm:%d %d %f|%d:%f,%f\n",iMP,iElm,dotProduct,i+1,edgeCenter[0],edgeCenter[1]);
                        }
                    if(dotProduct <=0){
                        edgeIndex = -1;
                        break;
                    }
                    if(currentDistSq < minDistSq){
                        edgeIndex = i+1;
                        minDistSq = currentDistSq;    
                    }
                }
                if(edgeIndex <0){
                    //we get to the final elm
                    MPs2Elm(iMP) = iElm;
                    MPsPosition(iMP) = MPnew;
                    break;
                }else{
                    //update the iELm and do the loop again
                    iElm = elm2ElmConn(iElm,edgeIndex);
                    //if(printVTP>=0)
                    //   Kokkos::atomic_increment(&count(iMP));
                }
            } 
        }
    });
    if(printVTP>=0){
        const int maxNum =5;
        IntView countNum("countNumCrossMPs",maxNum+1);
        Kokkos::parallel_for("countMPs",numMPs,KOKKOS_LAMBDA(const int iMP){
            for(int i=0; i<=maxNum; i++){
                if(i == count(iMP) || (i == maxNum && count(iMP)>= i)){
                    Kokkos::atomic_increment(&countNum(i));
                    break;
                }    
            }
        });
        IntView::HostMirror h_countNum = Kokkos::create_mirror_view(countNum);
        Kokkos::deep_copy(h_countNum, countNum);
        Kokkos::fence();
        for(int iCountNum = 0; iCountNum <= maxNum; iCountNum++){
            numMPs = h_countNum(iCountNum);
            printf("%d-%d:%d\n",iCountNum,printVTP,numMPs);
        }
    }
}

void MPM::CVTTrackingElmCenterBased(Vector2View dx, const int printVTP){
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

    //numMPs = 1; //XXX
    IntView count("countCrossMPs",numMPs);
    Vector2View elmCenter("elmentCenter",numElms);
    Kokkos::parallel_for("calcElmCenter",numElms,KOKKOS_LAMBDA(const int iElm){
        elmCenter(iElm) = calcElmCenter(iElm,elm2VtxConn,vtxCoords);
    });
    //printf("%d\n",numMPs);
    Kokkos::parallel_for("CVTElmCalc",numMPs,KOKKOS_LAMBDA(const int iMP){
        Vector2 MP = MPsPosition(iMP);
        if(isActive(iMP)){
            int iElm = MPs2Elm(iMP);
            Vector2 MPnew = MP + dx(iMP);
            //printf("%d:%f %f + %f %f = %f %f\n",iMP,MP[0],MP[1],dx(iMP)[0],dx(iMP)[1],MPnew[0],MPnew[1]);
            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                Vector2 delta = MPnew - elmCenter(iElm);
                double minDistSq = delta[0]*delta[0] + delta[1]*delta[1];
                int closestElm = -1;
                //go through all the connected elm, calc distance
                for(int i=1; i<=numVtx; i++){
                    int elmID = elm2ElmConn(iElm,i);
                    delta = MPnew - elmCenter(elmID);
                    double neighborDistSq = delta[0]*delta[0] + delta[1]*delta[1];
                    if(neighborDistSq < minDistSq){
                        closestElm = elmID;
                        minDistSq = neighborDistSq;
                    }
                }
                if(closestElm<0){
                    MPs2Elm(iMP) = iElm;
                    MPsPosition(iMP) = MPnew;
                    break;
                }else{
                    iElm = closestElm;
                    //if(printVTP>=0)
                    //   Kokkos::atomic_increment(&count(iMP));
                }
                //printf("ElmCenter: %d (%f,%f)\n",iElm,elmCenter(iElm)[0],elmCenter(iElm)[1]);
            } 
        }
    });
    if(printVTP>=0){
        const int maxNum =11;
        IntView countNum("countNumCrossMPs",maxNum+1);
        Kokkos::parallel_for("countMPs",numMPs,KOKKOS_LAMBDA(const int iMP){
            for(int i=0; i<=maxNum; i++){
                if(i == count(iMP) || (i == maxNum && count(iMP)>= i)){
                    Kokkos::atomic_increment(&countNum(i));
                    break;
                }    
            }
        });
        IntView::HostMirror h_countNum = Kokkos::create_mirror_view(countNum);
        Kokkos::deep_copy(h_countNum, countNum);
        Kokkos::fence();
        for(int iCountNum = 0; iCountNum <= maxNum; iCountNum++){
            numMPs = h_countNum(iCountNum);
            printf("%d-%d:%d\n",iCountNum,printVTP,numMPs);
        }
    }
   
}

void MPM::T2LTracking(Vector2View dx, const int printVTP){
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
   
    //numMPs = 100;//XXX
    Vector2View history("positionHistory",numMPs);
    Vector2View resultLeft("positionResult",numMPs);
    Vector2View resultRight("positionResult",numMPs); 
    IntView count("countCrossMPs",numMPs);
    Kokkos::parallel_for("T2LCalc",numMPs,KOKKOS_LAMBDA(const int iMP){
        Vector2 MP = MPsPosition(iMP);
        if(isActive(iMP)){
            int iElm = MPs2Elm(iMP);
            Vector2 MPnew = MP + dx(iMP);    
            
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
                        //go to the next elm
                        iElm = elm2ElmConn(iElm,i+1);
                        //if(printVTP>=0)
                        //    Kokkos::atomic_increment(&count(iMP));
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
                //if(printVTP>=0){ 
                //    Vector2 MParrow = MP + dx(iMP)*0.7;
                //    Vector2 shift = Vector2(-dx(iMP)[1],dx(iMP)[0])*0.1;
                //    Vector2 MPLeft = MParrow + shift;
                //    Vector2 MPRight = MParrow - shift;
                //   history(iMP) = MP;
                //    resultLeft(iMP) = MPLeft;
                //    resultRight(iMP) = MPRight;
                //}
                MPs2Elm(iMP) = iElm;
                MPsPosition(iMP) = MPnew;
                break;
            }
        }
        else{
            history(iMP) = MP;
            resultLeft(iMP) = MP;
            resultRight(iMP) = MP;
        }
    }); 
    if(printVTP>=0){
        //TODO: figure out the maxNum (parallel_reduce)
        const int maxNum = 5;

        IntView::HostMirror h_count = Kokkos::create_mirror_view(count);
        Vector2View::HostMirror h_history = Kokkos::create_mirror_view(history);
        Vector2View::HostMirror h_resultLeft = Kokkos::create_mirror_view(resultLeft);
        Vector2View::HostMirror h_resultRight = Kokkos::create_mirror_view(resultRight);
        Vector2View::HostMirror h_MPsPosition = Kokkos::create_mirror_view(MPsPosition);
        Kokkos::deep_copy(h_count,count);
        Kokkos::fence();
        IntView countNum("countNumCrossMPs",maxNum+1);
        Kokkos::parallel_for("countMPs",numMPs,KOKKOS_LAMBDA(const int iMP){
            for(int i=0; i<=maxNum; i++){
                if(i == count(iMP) || (i == maxNum && count(iMP)>= i)){
                    Kokkos::atomic_increment(&countNum(i));
                    break;
                }    
            }
        });
        Kokkos::deep_copy(h_history, history);
        Kokkos::deep_copy(h_resultLeft, resultLeft);
        Kokkos::deep_copy(h_resultRight, resultRight);
        Kokkos::deep_copy(h_MPsPosition, MPsPosition);
        IntView::HostMirror h_countNum = Kokkos::create_mirror_view(countNum); 
        Kokkos::deep_copy(h_countNum, countNum);
        Kokkos::fence();
        const int totalNumMPs = numMPs;
        for(int iCountNum = 0; iCountNum <= maxNum; iCountNum++){
            numMPs = h_countNum(iCountNum);
            printf("%d-%d:%d\n",iCountNum,printVTP,numMPs); 
//* printVTP file
            char* fileOutput = (char *)malloc(sizeof(char) * 256); 
            sprintf(fileOutput, "polyMpmTestVTPOutput_across%d-%d.vtp",iCountNum,printVTP);
            FILE * pFile = fopen(fileOutput,"w");
            free(fileOutput);   
            fprintf(pFile, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"%d\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",numMPs*4,numMPs*2); 
            for(int i=0; i<totalNumMPs; i++){
                if(h_count(i) == iCountNum || (iCountNum == maxNum && h_count(i) >= maxNum))
                //XXX: MPsPosition is the updated new position, h_history is the old position
                    fprintf(pFile,"          %f %f 0.0\n          %f %f 0.0\n          %f %f 0.0\n          %f %f 0.0\n",h_history(i)[0], h_history(i)[1], h_MPsPosition(i)[0], h_MPsPosition(i)[1], h_resultLeft(i)[0], h_resultLeft(i)[1], h_resultRight(i)[0], h_resultRight(i)[1]);
            }
            fprintf(pFile,"        </DataArray>\n      </Points>\n      <Lines>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"); 
            for(int i=0; i<numMPs*4; i+=4){
                // 01 213
                fprintf(pFile,"          %d %d\n          %d %d %d\n", i, i+1, i+2, i+1, i+3);
            }
            fprintf(pFile,"        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
            for(int i=0; i<numMPs*5; i+=5){
                fprintf(pFile,"          %d\n          %d\n",i+2,i+5);
            }
            fprintf(pFile,"        </DataArray>\n      </Lines>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
            fclose(pFile);
//===*/
        }
    }
}

} 
