#include <Kokkos_Core.hpp>
#include "pmt_utils.hpp"
#include "pmt_MPMesh.hpp"

namespace polyMpmTest{

void MPMesh::CVTTrackingEdgeCenterBased(Vector2View dx){
    auto mesh = getMesh();
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    auto numMPs = MPs->getCount();

    auto elm2VtxConn = mesh.getElm2VtxConn();
    auto vtx2ElmConn = mesh.getVtx2ElmConn();
    auto elm2ElmConn = mesh.getElm2ElmConn();
    const auto vtxCoords = mesh.getVtxCoords(); 
    auto mpPositions = MPs->getData<MPF_Cur_Pos_XYZ>();
    Kokkos::View<Vector2*[maxVtxsPerElm]> edgeCenters("EdgeCenters",numElms);
    
    auto calcCenter = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
        if(mask){
            int numVtx = elm2VtxConn(elm,0);
            int v[maxVtxsPerElm];
            for(int i=0; i< numVtx; i++)
                v[i] = elm2VtxConn(elm,i+1)-1;
            for(int i=0; i< numVtx; i++){
                Vector2 v_i = vtxCoords(v[i]);
                Vector2 v_ip1 = vtxCoords(v[(i+1)%numVtx]);
                edgeCenters(elm,i) = (v_ip1 + v_i)*0.5;
            }
        }
    };
    MPs->parallel_for(calcCenter,"calcEdgeCenter");
   
    auto CVTEdgeTracking = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
        Vector2 MP = Vector2(mpPositions(mp,0),mpPositions(mp,1));//XXX:the input is XYZ, but we only support 2d vector
        if(mask){
            Vector2 MPnew = MP + dx(mp);
            int iElm = elm;
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
                    double dotProduct = dx(mp).dot(delta);
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
                    //MPs2Elm(mp) = elm; TODO:find a to set the new MP position
                    mpPositions(mp,0) = MPnew[0];
                    mpPositions(mp,1) = MPnew[1];
                    mpPositions(mp,2) = 0.0; //XXX:we only have 2d vector
                    break;
                }else{
                    //update the iELm and do the loop again
                    iElm = elm2ElmConn(iElm,edgeIndex);
                }
            } 
        }
    };
    MPs->parallel_for(CVTEdgeTracking,"CVTTrackingEdgeCenterBasedCalc");
}


void MPMesh::CVTTrackingElmCenterBased(Vector2View dx){
    int numVtxs = mesh_.getNumVertices();
    int numElms = mesh_.getNumElements();
    
    const auto vtxCoords = mesh_.getVtxCoords(); 
    auto elm2VtxConn = mesh_.getElm2VtxConn();
    auto vtx2ElmConn = mesh_.getVtx2ElmConn();
    auto elm2ElmConn = mesh_.getElm2ElmConn();

    auto mpPositions = MPs->getData<MPF_Cur_Pos_XYZ>();
    //auto MPs2Elm = materialPoints2Elm_;

    Vector2View elmCenter("elmentCenter",numElms);
    auto calcCenter = PS_LAMBDA(const int& elm, const int& mp, const int&mask){
//        elmCenter(iElm) = calcElmCenter(iElm,elm2VtxConn,vtxCoords);
        int numVtx = elm2VtxConn(elm,0);
        double sum_x = 0.0, sum_y = 0.0;
        for(int i=1; i<= numVtx; i++){
            sum_x += vtxCoords(elm2VtxConn(elm,i)-1)[0];
            sum_y += vtxCoords(elm2VtxConn(elm,i)-1)[1];
        }
    };
    MPs->parallel_for(calcCenter,"calcElmCenter");

    auto CVTElmCalc = PS_LAMBDA(const int& elm, const int& mp, const int&mask){
        Vector2 MP = Vector2(mpPositions(mp,0),mpPositions(mp,1));//XXX:the input is XYZ, but we only support 2d vector
        if(mask){
            int iElm = elm;
            Vector2 MPnew = MP + dx(mp);
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
                    //MPs2Elm(mp) = iElm;
                    mpPositions(mp,0) = MPnew[0];
                    mpPositions(mp,1) = MPnew[1];
                    mpPositions(mp,2) = 0.0; //XXX:we only have 2d vector
                    break;
                }else{
                    iElm = closestElm;
                }
            } 
        }
    };
    MPs->parallel_for(CVTElmCalc,"CVTTrackingElmCenterBasedCalc");
}

void MPMesh::T2LTracking(Vector2View dx){
    int numVtxs = mesh_.getNumVertices();
    int numElms = mesh_.getNumElements();
    
    const auto vtxCoords = mesh_.getVtxCoords(); 
    auto elm2VtxConn = mesh_.getElm2VtxConn();
    auto vtx2ElmConn = mesh_.getVtx2ElmConn();
    auto elm2ElmConn = mesh_.getElm2ElmConn();

    auto mpPositions = MPs->getData<MPF_Cur_Pos_XYZ>();
    //auto MPs2Elm = materialPoints2Elm_;
   
    auto T2LCalc = PS_LAMBDA(const int& elm, const int& mp, const int&mask){
        Vector2 MP = Vector2(mpPositions(mp,0),mpPositions(mp,1));//XXX:the input is XYZ, but we only support 2d vector
        if(mask){
            int iElm = elm;
            Vector2 MPnew = MP + dx(mp);    
            
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
                    pdx[i] = (v_i - MP).cross(dx(mp));
                }
                
                for(int i=0; i<numVtx; i++){
                    int ip1 = (i+1)%numVtx;
                    //pdx*pdx<0 and
                    if(pdx[i]*pdx[ip1] <0 && e[i].cross(MPnew-vtxCoords(v[i]))<0){
                        //go to the next elm
                        iElm = elm2ElmConn(iElm,i+1);
                        goToNeighbour = true;
                        if(iElm <0){
                            //mask = false;
                            goToNeighbour = false;
                        }
                    }
                }
                //if goes to the other 
                if(goToNeighbour)
                    continue; 
                //otherwise we do the update and end the loop
                //MPs2Elm(iMP) = iElm; TODO
                mpPositions(mp,0) = MPnew[0];
                mpPositions(mp,1) = MPnew[1];
                mpPositions(mp,2) = 0.0; //XXX:we only have 2d vector
                break;
            }
        }
    }; 
    MPs->parallel_for(T2LCalc,"T2lTrackingCalc");
}

} 
