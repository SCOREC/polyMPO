#include <Kokkos_Random.hpp>
#include "testUtils.hpp"

using namespace polyMpmTest;

void interpolateWachspress(MPMesh& mpMesh){
    auto mesh = mpMesh.getMesh();
    auto vtxCoords = mesh.getVtxCoords();
    auto elm2VtxConn = mesh.getElm2VtxConn();

    auto MPs = mpMesh.MPs;
    auto MPsPosition = MPs->getPositions();
    auto eval = PS_LAMBDA(const int& elm, const int& mp, const int mask){
        if (mask) {
            Vector2 v[maxVtxsPerElm+1] = {vtxCoords(elm2VtxConn(elm,1))};
            initArray(v,maxVtxsPerElm+1,vtxCoords(elm2VtxConn(elm,1)));
            int numVtx = elm2VtxConn(elm,0);
            for(int i = 1; i<=numVtx; i++){
                v[i-1] = vtxCoords(elm2VtxConn(elm,i)-1);
            }
            v[numVtx] = vtxCoords(elm2VtxConn(elm,1)-1);
            double basisByArea[maxVtxsPerElm] = {0.0};
            initArray(basisByArea,maxVtxsPerElm,0.0);
            Vector2 gradBasisByArea[maxVtxsPerElm];
            Vector2 position(MPsPosition(mp,0),MPsPosition(mp,1));
            getBasisAndGradByAreaGblForm(position, numVtx, v, basisByArea, gradBasisByArea);
        
            Vector2 wp_coord(0.0,0.0);
            double wp_grad = 0.0;
            for(int i=0; i<= numVtx; i++){
                wp_coord = wp_coord + v[i]*basisByArea[i];
                wp_grad = wp_grad + gradBasisByArea[i].dot(v[i]);
            }
        }        
    };
    MPs->parallel_for(eval, "getBasisAndGradByAreaGblForm");
}
