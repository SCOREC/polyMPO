#include <Kokkos_Random.hpp>
#include "testUtils.hpp"

using namespace polyMPO;

void interpolateWachspress(MPMesh& mpMesh){
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getVtxCoords();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();

    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    auto eval = PS_LAMBDA(const int& elm, const int& mp, const int mask){
        if (mask) {
            //convert the double[] to Vec2d 
            Vec2d v[maxVtxsPerElm+1];
            Vec3d v3d[maxVtxsPerElm+1];
            initArray(v,maxVtxsPerElm+1,Vec2d());
            int numVtx = elm2VtxConn(elm,0);
            for(int i = 1; i<=numVtx; i++){
                v[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
                v[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
                v3d[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
                v3d[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
                v3d[i-1][2] = vtxCoords(elm2VtxConn(elm,i)-1,2);
            }
            v[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
            v[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
            v3d[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
            v3d[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
            v3d[numVtx][2] = vtxCoords(elm2VtxConn(elm,1)-1,2);
            double basisByArea[maxVtxsPerElm] = {0.0};
            initArray(basisByArea,maxVtxsPerElm,0.0);
            Vec2d gradBasisByArea[maxVtxsPerElm];
            Vec2d position(MPsPosition(mp,0),MPsPosition(mp,1));
            getBasisAndGradByAreaGblForm(position, numVtx, v, basisByArea, gradBasisByArea);

            double basisByArea3d[maxVtxsPerElm] = {0.0};
            initArray(basisByArea3d,maxVtxsPerElm,0.0);
            double radius = 1.0;
            Vec3d position3d(MPsPosition(mp,0),MPsPosition(mp,1),MPsPosition(mp,2));
            getBasisByAreaGblFormSpherical(position3d, numVtx, v3d, radius, basisByArea3d);
        
            Vec2d wp_coord(0.0,0.0);
            double wp_grad = 0.0;
            for(int i=0; i<= numVtx; i++){
                wp_coord = wp_coord + v[i]*basisByArea[i];
                wp_grad = wp_grad + gradBasisByArea[i].dot(v[i]);
            }
        }        
    };
    p_MPs->parallel_for(eval, "getBasisAndGradByAreaGblForm");
}
