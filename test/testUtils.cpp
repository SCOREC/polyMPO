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
                //printf("%d:(%f,%f,%f)\n",i-1,v3d[i-1][0],v3d[i-1][1],v3d[i-1][2]);
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

void interpolateWachspress3DTest(MPMesh& mpMesh){
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();

    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();
    auto eval = PS_LAMBDA(const int& elm, const int& mp, const int mask){
        if (mask && mp == 1) {
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
           
            // testing
            /*const int numVtx_test = 4;
            Vec3d position_test(3.0, 2.0, 3.0); // MP
            Vec3d v_test[numVtx_test];
            double basisByArea_test[numVtx_test+1];
            Vec3d gradBasisByArea_test[numVtx_test];
            v_test[0] = Vec3d(3.0, 2.0, 5.0);
            v_test[1] = Vec3d(3.0, 3.0, 2.0);
            v_test[2] = Vec3d(2.0, 1.0, 2.0);
            v_test[3] = Vec3d(5.0, 1.0, 2.0);
            v_test[4] = v_test[0];
            getBasisAndGradByAreaGblForm3d(position_test, numVtx_test, v_test, basisByArea_test, gradBasisByArea_test);

            for (int i = 0; i <numVtx_test; i++) {
                printf("%d: %f %f \n", i, gradBasisByArea_test[i][0],  gradBasisByArea_test[i][1]);
            }
            */
            // testing ends
            double rt = 1.41421356237; // sqrt 2
            double rth = 1.73205080757;// sqrt 3

            // rotation matrix
            Vec3d r[3] = {Vec3d(1.0/rt,-rth/(2*rt),1/(2*rt)),
                          Vec3d(1.0/rt, rth/(2*rt),-1.0/(2*rt)),
                          Vec3d(0, 1.0/2.0, rth/2.0)};
            // r inverse
            Vec3d ri[3] = {Vec3d(rt - 1.0/rt,1/(rt),0),
                          Vec3d(-rt/rth + 1.0/(2.0*rt*rth), rt/rth - 1.0/(2.0*rt*rth),1.0/2.0),
                          Vec3d(1.0/(2.0*rt), -1.0/(2.0*rt), rth/2.0)};
    
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
                //double fi = af * v[i][0] + bf * v[i][1] + cf * v[i][2] + kf;
                double v1 = v[i].dot(ri[0]);
                double v2 = v[i].dot(ri[1]);
                double v3 = v[i].dot(ri[2]);
                double fi = af * v1 + bf * v2 + cf * v3 + kf;
                wp_grad = wp_grad + gradBasisByArea[i] * fi;
                /*printf("i: %d, gradBasis: %.16e %.16e %.16e fi: %.16e \nx: %.16e y: %.16e z: %.16e \n", i,
                                                                              gradBasisByArea[i][0],
                                                                              gradBasisByArea[i][1],
                                                                              gradBasisByArea[i][2],
                                                                              fi,
                                                                              v1,
                                                                              v2,
                                                                              v3); 
                */
            }
              
            double gxt = af; // gx tilde
            double gyt = bf; // gy tilde

            Vec3d wp_grad2;
            wp_grad2[0] = gxt * r[0][0] + gyt * r[0][1];
            wp_grad2[1] = gxt * r[1][0] + gyt * r[1][1];
            wp_grad2[2] = gxt * r[2][0] + gyt * r[2][1];
           
            printf("WP gradient:(%.16e %.16e %.16e)\nexpected gradient:(%.16e %.16e %.16e)\n",
                                              wp_grad[0],
                                              wp_grad[1],
                                              wp_grad[2],
                                              wp_grad2[0],
                                              wp_grad2[1],
                                              wp_grad2[2]);
            
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

void printVTP(MPMesh& mpMesh){
    auto p_mesh = mpMesh.p_mesh;
    auto vtxCoords = p_mesh->getMeshField<polyMPO::MeshF_VtxCoords>();
    auto elm2VtxConn = p_mesh->getElm2VtxConn();

    auto p_MPs = mpMesh.p_MPs;
    auto MPsPosition = p_MPs->getPositions();

    DoubleVec3dView::HostMirror h_vtxCoords = Kokkos::create_mirror_view(vtxCoords);
    IntVtx2ElmView::HostMirror h_elm2VtxConn = Kokkos::create_mirror_view(elm2VtxConn);
    const int nCells = p_mesh->getNumElements();
    const int nVertices = p_mesh->getNumVertices();
    Kokkos::deep_copy(h_vtxCoords,vtxCoords);
    Kokkos::deep_copy(h_elm2VtxConn,elm2VtxConn);
    printf("<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n  <PolyData>\n    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%d\">\n      <Points>\n        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n",nVertices,nCells);
    for(int i=0; i<nVertices; i++){
        printf("          %f %f %f\n",h_vtxCoords(i,0),h_vtxCoords(i,1),h_vtxCoords(i,2));
    }
    printf("        </DataArray>\n      </Points>\n      <Polys>\n        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for(int i=0; i<nCells; i++){
        printf("          ");
        for(int j=0; j< h_elm2VtxConn(i,0); j++){
            printf("%d ", h_elm2VtxConn(i,j+1)-1);
        } 
        printf("\n");
    }
    printf("        </DataArray>\n        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
    
    int count = 0;
    for(int i=0;i<nCells; i++){
        count += h_elm2VtxConn(i,0);
        printf("          %d\n",count);
    }
    printf("        </DataArray>\n      </Polys>\n    </Piece>\n  </PolyData>\n</VTKFile>\n");
}

void readMPAS(){}
