#include <Kokkos_Core.hpp>
#include "pmpo_utils.hpp"
#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"

namespace polyMPO{

void printVTP_mesh(MPMesh& mpMesh, int printVTPIndex=-1);

void MPMesh::CVTTrackingEdgeCenterBased(Vec2dView dx){
    int numElms = p_mesh->getNumElements();

    auto elm2VtxConn = p_mesh->getElm2VtxConn();
    auto elm2ElmConn = p_mesh->getElm2ElmConn();
    auto MPs2Elm = p_MPs->getData<MPF_Tgt_Elm_ID>();
    const auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>(); 
    auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
    Kokkos::View<Vec2d*[maxVtxsPerElm]> edgeCenters("EdgeCenters",numElms);
   
    Kokkos::parallel_for("calcEdgeCenter", numElms, KOKKOS_LAMBDA(const int elm){  
        int numVtx = elm2VtxConn(elm,0);
        int v[maxVtxsPerElm];
        for(int i=0; i< numVtx; i++)
            v[i] = elm2VtxConn(elm,i+1)-1;
        for(int i=0; i< numVtx; i++){
            int idx_ip1 = (i+1)%numVtx;
            Vec2d v_i(vtxCoords(v[i],0),vtxCoords(v[i],1));
            Vec2d v_ip1(vtxCoords(v[idx_ip1],0),vtxCoords(v[idx_ip1],1));
            edgeCenters(elm,i) = (v_ip1 + v_i)*0.5;
        }
    });
   
    auto CVTEdgeTracking = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
        Vec2d MP(mpPositions(mp,0),mpPositions(mp,1));//XXX:the input is XYZ, but we only support 2d vector
        if(mask){
            Vec2d MPnew = MP + dx(mp);
            int iElm = elm;
            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                //calc dist square from each edge center to MPnew
                //calc dot products to check inside or not
                int edgeIndex = -1;
                double minDistSq = DBL_MAX;
                for(int i=0; i< numVtx; i++){
                    Vec2d edgeCenter = edgeCenters(iElm,i);
                    Vec2d delta = MPnew - edgeCenter;
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
                    MPs2Elm(mp) = iElm; 
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
    p_MPs->parallel_for(CVTEdgeTracking,"CVTTrackingEdgeCenterBasedCalc");
}


void MPMesh::CVTTrackingElmCenterBased(const int printVTPIndex){
    int numVtxs = p_mesh->getNumVertices();
    int numElms = p_mesh->getNumElements();
    auto numMPs = p_MPs->getCount();

    const auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>(); 
    auto elm2VtxConn = p_mesh->getElm2VtxConn();
    auto elm2ElmConn = p_mesh->getElm2ElmConn();

    auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
    auto mpTgtPos = p_MPs->getData<MPF_Tgt_Pos_XYZ>();
    auto MPs2Elm = p_MPs->getData<MPF_Tgt_Elm_ID>();;
    
    if(printVTPIndex>=0) {
      printVTP_mesh(printVTPIndex);
    }

    Vec3dView elmCenter("elementCenter",numElms);
    Kokkos::parallel_for("calcElementCenter", numElms, KOKKOS_LAMBDA(const int elm){  
        int numVtx = elm2VtxConn(elm,0);
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
        for(int i=1; i<= numVtx; i++){
            sum_x += vtxCoords(elm2VtxConn(elm,i)-1,0);
            sum_y += vtxCoords(elm2VtxConn(elm,i)-1,1);
            sum_z += vtxCoords(elm2VtxConn(elm,i)-1,2);
        }
        elmCenter(elm)[0] = sum_x/numVtx;
        elmCenter(elm)[1] = sum_y/numVtx;
        elmCenter(elm)[2] = sum_z/numVtx;
    });

    Vec3dView history("positionHistory",numMPs);
    Vec3dView resultLeft("positionResult",numMPs);
    Vec3dView resultRight("positionResult",numMPs);
    Vec3dView mpTgtPosArray("positionTarget",numMPs);

    auto CVTElmCalc = PS_LAMBDA(const int& elm, const int& mp, const int&mask){
        Vec3d MP(mpPositions(mp,0),mpPositions(mp,1),mpPositions(mp,2));
        if(mask){
            int iElm = elm;
            Vec3d MPnew(mpTgtPos(mp,0),mpTgtPos(mp,1),mpTgtPos(mp,2));
            Vec3d dx = MPnew-MP;
            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                Vec3d delta = MPnew - elmCenter(iElm);
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
                    MPs2Elm(mp) = iElm;
                    break;
                }else{
                    iElm = closestElm;
                }
            }
            if(printVTPIndex>=0){ 
                double d1 = dx[0];
                double d2 = dx[2];
                double d3 = dx[3];
                double m1 = MP[0];
                double m2 = MP[1];
                double m3 = MP[2];
                Vec3d MParrow = MP + dx*0.7;
                Vec3d r = MPnew * (1.0/MPnew.magnitude());
                Vec3d shift = dx.cross(r) * ((1.0-0.7)*dx.magnitude()/(dx.cross(r)).magnitude());
                Vec3d MPLeft = MParrow + shift;
                Vec3d MPRight = MParrow - shift;
                history(mp) = MP;
                resultLeft(mp) = MPLeft;
                resultRight(mp) = MPRight;
                mpTgtPosArray(mp) = MPnew;
            }
        }
    };
    p_MPs->parallel_for(CVTElmCalc,"CVTTrackingElmCenterBasedCalc");

    if(printVTPIndex>=0){
        Vec3dView::HostMirror h_history = Kokkos::create_mirror_view(history);
        Vec3dView::HostMirror h_resultLeft = Kokkos::create_mirror_view(resultLeft);
        Vec3dView::HostMirror h_resultRight = Kokkos::create_mirror_view(resultRight);
        Vec3dView::HostMirror h_mpTgtPos = Kokkos::create_mirror_view(mpTgtPosArray);

        Kokkos::deep_copy(h_history, history);
        Kokkos::deep_copy(h_resultLeft, resultLeft);
        Kokkos::deep_copy(h_resultRight, resultRight);
        Kokkos::deep_copy(h_mpTgtPos, mpTgtPosArray);

        // printVTP file
        char* fileOutput = (char *)malloc(sizeof(char) * 256); 
        sprintf(fileOutput, "polyMPOCVTTrackingElmCenter_MPtracks_%d.vtp", printVTPIndex);
        FILE * pFile = fopen(fileOutput,"w");
        free(fileOutput);   
        fprintf(pFile, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"%d\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",numMPs*4,numMPs*2); 
        for(int i=0; i<numMPs; i++){
            fprintf(pFile,"          %f %f %f\n          %f %f %f\n          %f %f %f\n          %f %f %f\n",
                          h_history(i)[0],h_history(i)[1],h_history(i)[2],
                          h_mpTgtPos(i)[0],h_mpTgtPos(i)[1],h_mpTgtPos(i)[2],
                          h_resultLeft(i)[0],h_resultLeft(i)[1],h_resultLeft(i)[2],
                          h_resultRight(i)[0],h_resultRight(i)[1],h_resultRight(i)[2]);
        }
        fprintf(pFile,"        </DataArray>\n      </Points>\n      <Lines>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n"); 
        for(int i=0; i<numMPs*4; i+=4){
             fprintf(pFile,"          %d %d\n          %d %d %d\n",i,i+1,i+2,i+1,i+3);
        }
        fprintf(pFile,"        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        for(int i=0; i<numMPs*5; i+=5){
            fprintf(pFile,"          %d\n          %d\n",i+2,i+5);
        }
        fprintf(pFile,"        </DataArray>\n      </Lines>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
        fclose(pFile);
    }
}

void MPMesh::T2LTracking(Vec2dView dx){
    const auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>(); 
    auto elm2VtxConn = p_mesh->getElm2VtxConn();
    auto elm2ElmConn = p_mesh->getElm2ElmConn();

    auto mpPositions = p_MPs->getData<MPF_Cur_Pos_XYZ>();
    auto MPs2Elm = p_MPs->getData<MPF_Tgt_Elm_ID>();
    auto mpStatus = p_MPs->getData<MPF_Status>();
   
    auto T2LCalc = PS_LAMBDA(const int& elm, const int& mp, const int&mask){
        Vec2d MP(mpPositions(mp,0),mpPositions(mp,1));//XXX:the input is XYZ, but we only support 2d vector
        if(mask){
            int iElm = elm;
            Vec2d MPnew = MP + dx(mp);    
            
            while(true){
                int numVtx = elm2VtxConn(iElm,0);
                bool goToNeighbour = false;
                //seperate the elm2Vtx
                int v[maxVtxsPerElm];
                for(int i=0; i< numVtx; i++)
                    v[i] = elm2VtxConn(iElm,i+1)-1;
                //get edges and perpendiculardx
                Vec2d e[maxVtxsPerElm];
                double pdx[maxVtxsPerElm];                    
                for(int i=0; i< numVtx; i++){
                    int idx_ip1 = (i+1)%numVtx;
                    Vec2d v_i(vtxCoords(v[i],0),vtxCoords(v[i],1));
                    Vec2d v_ip1(vtxCoords(v[idx_ip1],0),vtxCoords(v[idx_ip1],1));
                    e[i] = v_ip1 - v_i;
                    pdx[i] = (v_i - MP).cross(dx(mp));
                }
                
                for(int i=0; i<numVtx; i++){
                    int ip1 = (i+1)%numVtx;
                    //pdx*pdx<0 and edge is acrossed 
                    if(pdx[i]*pdx[ip1] <0 && e[i].cross(Vec2d(MPnew[0]-vtxCoords(v[i],0),
                                                              MPnew[1]-vtxCoords(v[i],1)))<0){
                        //go to the next elm
                        iElm = elm2ElmConn(iElm,i+1);
                        goToNeighbour = true;
                        if(iElm <0){
                            mpStatus(mp) = 0;                  
                            MPs2Elm(mp) = -1;
                            goToNeighbour = false;
                        }
                    }
                }
                //if goes to the other 
                if(goToNeighbour)
                    continue; 
                //otherwise we do the update and end the loop
                MPs2Elm(mp) = iElm;
                mpPositions(mp,0) = MPnew[0];
                mpPositions(mp,1) = MPnew[1];
                mpPositions(mp,2) = 0.0; //XXX:we only have 2d vector
                break;
            }
        }
    }; 
    p_MPs->parallel_for(T2LCalc,"T2lTrackingCalc");
}

void MPMesh::reconstructSlices() {
    // if (!reconsOption) return;
    p_MPs->calcBasis();
    for (auto const& [index, reconstruct] : reconstructSlice) {
        if (reconstruct) reconstruct();
    }
    reconstructSlice.clear();
}

void MPMesh::push(){
  p_mesh->computeRotLatLonIncr();
  sphericalInterpolation<MeshF_RotLatLonIncr>(*this);
  p_MPs->updateRotLatLonAndXYZ2Tgt(p_mesh->getSphereRadius()); // set Tgt_XYZ

  CVTTrackingElmCenterBased(); // move to Tgt_XYZ
  reconstructSlices();

  p_MPs->updateMPSlice<MPF_Cur_Pos_XYZ, MPF_Tgt_Pos_XYZ>(); // Tgt_XYZ becomes Cur_XYZ
  p_MPs->updateMPSlice<MPF_Cur_Pos_Rot_Lat_Lon, MPF_Tgt_Pos_Rot_Lat_Lon>(); // Tgt becomes Cur
  p_MPs->rebuild(); //rebuild pumi-pic
  p_MPs->updateMPElmID(); //update mpElm IDs slices
}

void MPMesh::printVTP_mesh(int printVTPIndex){
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();

    auto MPsPosition = p_MPs->getPositions();

    char* fileOutput = (char *)malloc(sizeof(char) * 256); 
    sprintf(fileOutput,"polyMPO_MPMesh_mesh_%d.vtp", printVTPIndex);
    FILE * pFile = fopen(fileOutput,"w");
    free(fileOutput);

    auto h_vtxCoords = Kokkos::create_mirror_view(vtxCoords);
    IntVtx2ElmView::HostMirror h_elm2VtxConn = Kokkos::create_mirror_view(elm2VtxConn);
    const int nCells = p_mesh->getNumElements();
    const int nVertices = p_mesh->getNumVertices();
    Kokkos::deep_copy(h_vtxCoords,vtxCoords);
    Kokkos::deep_copy(h_elm2VtxConn,elm2VtxConn);
    fprintf(pFile, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%d\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",nVertices,nCells);
    for(int i=0; i<nVertices; i++){
        fprintf(pFile, "          %f %f %f\n",h_vtxCoords(i,0),h_vtxCoords(i,1),h_vtxCoords(i,2));
    }
    fprintf(pFile, "        </DataArray>\n      </Points>\n      <Polys>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int i=0; i<nCells; i++){
        fprintf(pFile, "          ");
        for(int j=0; j< h_elm2VtxConn(i,0); j++){
            fprintf(pFile, "%d ", h_elm2VtxConn(i,j+1)-1);
        } 
        fprintf(pFile, "\n");
    }
    fprintf(pFile, "        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    
    int count = 0;
    for(int i=0;i<nCells; i++){
        count += h_elm2VtxConn(i,0);
        fprintf(pFile, "          %d\n",count);
    }
    fprintf(pFile, "        </DataArray>\n      </Polys>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
    fclose(pFile);
}

}
