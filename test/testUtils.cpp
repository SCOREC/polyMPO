#include <Kokkos_Random.hpp>
#include "testUtils.hpp"

using namespace polyMPO;

void interpolateWachspress2DTest(MPMesh& mpMesh){
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();

    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    auto eval = PS_LAMBDA(const int& elm, const int& mp, const int mask){
        if (mask) {
            //convert the double[] to Vec2d 
            Vec2d v[maxVtxsPerElm+1];
            initArray(v,maxVtxsPerElm+1,Vec2d());
            int numVtx = elm2VtxConn(elm,0);
            for(int i = 1; i<=numVtx; i++){
                v[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
                v[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
            }
            v[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
            v[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
            double basisByArea[maxVtxsPerElm];
            initArray(basisByArea,maxVtxsPerElm,0.0);
            double basisByArea2[maxVtxsPerElm];
            initArray(basisByArea2,maxVtxsPerElm,0.0);
            Vec2d gradBasisByArea[maxVtxsPerElm];
            Vec2d position(MPsPosition(mp,0),MPsPosition(mp,1));
            getBasisAndGradByAreaGblForm2d(position, numVtx, v, basisByArea, gradBasisByArea);
            getBasisByAreaGblForm(position, numVtx, v, basisByArea2);

            double af = 10.1;
            double bf = 1.34;
            double kf = 3.55;
            Vec2d wp_coord(0.0,0.0);
            Vec2d wp_coord2(0.0,0.0);
            Vec2d wp_grad(0.0,0.0);
            for(int i=0; i< numVtx; i++){
                wp_coord = wp_coord + v[i]*basisByArea[i];
                wp_coord2 = wp_coord2 + v[i]*basisByArea2[i];
                double fi = af * v[i][0] + bf * v[i][1] + kf;
                wp_grad = wp_grad + gradBasisByArea[i] * fi;
            }
            assert(abs(wp_coord[0] - MPsPosition(mp,0)) < TEST_EPSILON);
            assert(abs(wp_coord[1] - MPsPosition(mp,1)) < TEST_EPSILON);
            assert(abs(wp_coord2[0] - MPsPosition(mp,0)) < TEST_EPSILON);
            assert(abs(wp_coord2[1] - MPsPosition(mp,1)) < TEST_EPSILON);
            assert(abs(wp_grad[0] - af) < TEST_EPSILON);
            assert(abs(wp_grad[1] - bf) < TEST_EPSILON);
        }        
    };
    p_MPs->parallel_for(eval, "interpolateWachspress2DTest");
}

void interpolateWachspressSphericalTest(MPMesh& mpMesh){
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();

    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    double radius = p_mesh->getSphereRadius();
    PMT_ALWAYS_ASSERT(radius >0);
    auto eval = PS_LAMBDA(const int& elm, const int& mp, const int mask){
        if (mask) {
            Vec3d position3d(MPsPosition(mp,0),MPsPosition(mp,1),MPsPosition(mp,2));
            Vec3d v3d[maxVtxsPerElm+1];
            int numVtx = elm2VtxConn(elm,0);
            for(int i = 1; i<=numVtx; i++){
                v3d[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
                v3d[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
                v3d[i-1][2] = vtxCoords(elm2VtxConn(elm,i)-1,2);
            }
            v3d[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
            v3d[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
            v3d[numVtx][2] = vtxCoords(elm2VtxConn(elm,1)-1,2);

            double basisByArea3d[maxVtxsPerElm] = {0.0};
            initArray(basisByArea3d,maxVtxsPerElm,0.0);
            getBasisByAreaGblFormSpherical(position3d, numVtx, v3d, radius, basisByArea3d);

            double basisByArea3d2[maxVtxsPerElm] = {0.0};
            initArray(basisByArea3d2,maxVtxsPerElm,0.0);
            getBasisByAreaGblFormSpherical2(position3d, numVtx, v3d, radius, basisByArea3d2);
       
            double basisByArea3d3[maxVtxsPerElm] = {0.0};
            initArray(basisByArea3d3,maxVtxsPerElm,0.0);
            getBasisByAreaGblFormSpherical3(position3d, numVtx, v3d, basisByArea3d3);
         
            Vec3d wp_coord(0.0,0.0,0.0);
            Vec3d wp_coord2(0.0,0.0,0.0);
            for(int i=0; i<numVtx; i++){
                wp_coord = wp_coord + v3d[i]*basisByArea3d[i];
                wp_coord2 = wp_coord + v3d[i]*basisByArea3d2[i];
            }
        }        
    };
    p_MPs->parallel_for(eval, "interpolateWachspressSphericalTest");
}

void interpolateWachspress3DTest(MPMesh& mpMesh, const int testMeshOption){
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();

    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    auto eval = PS_LAMBDA(const int& elm, const int& mp, const int mask){
        if (mask) {
            //convert the double[] to Vec3d 
            Vec3d v[maxVtxsPerElm+1];
            initArray(v,maxVtxsPerElm+1,Vec3d());
            int numVtx = elm2VtxConn(elm,0);
            for(int i = 1; i<=numVtx; i++){
                v[i-1][0] = vtxCoords(elm2VtxConn(elm,i)-1,0);
                v[i-1][1] = vtxCoords(elm2VtxConn(elm,i)-1,1);
                v[i-1][2] = vtxCoords(elm2VtxConn(elm,i)-1,2);
            }
            v[numVtx][0] = vtxCoords(elm2VtxConn(elm,1)-1,0);
            v[numVtx][1] = vtxCoords(elm2VtxConn(elm,1)-1,1);
            v[numVtx][2] = vtxCoords(elm2VtxConn(elm,1)-1,2);
            double basisByArea[maxVtxsPerElm];
            initArray(basisByArea,maxVtxsPerElm,0.0);
            double basisByArea2[maxVtxsPerElm];
            initArray(basisByArea2,maxVtxsPerElm,0.0);
            Vec3d gradBasisByArea[maxVtxsPerElm];
            Vec3d position(MPsPosition(mp,0),MPsPosition(mp,1),MPsPosition(mp,2));
            getBasisAndGradByAreaGblForm3d(position, numVtx, v, basisByArea, gradBasisByArea);
            getBasisByAreaGblForm3d(position, numVtx, v, basisByArea2);
          
            // rotation matrix 
            Vec3d r[3] = {Vec3d(1.0,0.0,0.0),
                          Vec3d(0.0,1.0,0.0),
                          Vec3d(0.0,0.0,1.0)};
            // r inverse
            Vec3d ri[3] = {Vec3d(1.0,0.0,0.0),
                           Vec3d(0.0,1.0,0.0),
                           Vec3d(0.0,0.0,1.0)};
            
            if (testMeshOption == 2) {
                // rotated 30 degrees around x-axis, then 45 degrees around z-axis
                double rt = 1.41421356237; // sqrt 2
                double rth = 1.73205080757;// sqrt 3
               
                // rotation matrix
                r[0] = Vec3d(1.0/rt,-rth/(2*rt),1/(2*rt));
                r[1] = Vec3d(1.0/rt, rth/(2*rt),-1.0/(2*rt));
                r[2] = Vec3d(0,1.0/2.0,rth/2.0);
                
                // r inverse
                ri[0] = Vec3d(rt - 1.0/rt,1.0/(rt),0.0);
                ri[1] = Vec3d(-rt/rth + 1.0/(2.0*rt*rth),rt/rth - 1.0/(2.0*rt*rth),1.0/2.0);
                ri[2] = Vec3d(1.0/(2.0*rt),-1.0/(2.0*rt),rth/2.0);
            }
            double af = 10.1;
            double bf = 1.34;
            double cf = 100.45;
            double kf = 3.55;
            Vec3d wp_coord(0.0,0.0,0.0);
            Vec3d wp_coord2(0.0,0.0,0.0);
            Vec3d wp_grad(0.0,0.0,0.0);
            for(int i=0; i<numVtx; i++){
                wp_coord = wp_coord + v[i]*basisByArea[i];
                wp_coord2 = wp_coord2 + v[i]*basisByArea2[i];
                double v1 = v[i].dot(ri[0]);
                double v2 = v[i].dot(ri[1]);
                double v3 = v[i].dot(ri[2]);
                double fi = af * v1 + bf * v2 + cf * v3 + kf;
                wp_grad = wp_grad + gradBasisByArea[i] * fi;
            }
              
            double gxt = af; // gx tilde
            double gyt = bf; // gy tilde
            Vec3d wp_grad2;
            wp_grad2[0] = gxt * r[0][0] + gyt * r[0][1];
            wp_grad2[1] = gxt * r[1][0] + gyt * r[1][1];
            wp_grad2[2] = gxt * r[2][0] + gyt * r[2][1];
            
            assert(abs(wp_coord[0] - MPsPosition(mp,0)) < TEST_EPSILON);
            assert(abs(wp_coord[1] - MPsPosition(mp,1)) < TEST_EPSILON);
            assert(abs(wp_coord[2] - MPsPosition(mp,2)) < TEST_EPSILON);
            assert(abs(wp_coord2[0] - MPsPosition(mp,0)) < TEST_EPSILON);
            assert(abs(wp_coord2[1] - MPsPosition(mp,1)) < TEST_EPSILON);
            assert(abs(wp_coord2[2] - MPsPosition(mp,2)) < TEST_EPSILON);
            assert(abs(wp_grad[0] - wp_grad2[0]) < TEST_EPSILON);
            assert(abs(wp_grad[1] - wp_grad2[1]) < TEST_EPSILON);
            assert(abs(wp_grad[2] - wp_grad2[2]) < TEST_EPSILON);
        }        
    };
    p_MPs->parallel_for(eval, "interpolateWachspress3DTest");
}

