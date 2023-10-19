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
            double basisByArea[maxVtxsPerElm] = {0.0};
            initArray(basisByArea,maxVtxsPerElm,0.0);
            Vec2d gradBasisByArea[maxVtxsPerElm];
            Vec2d position(MPsPosition(mp,0),MPsPosition(mp,1));
            getBasisAndGradByAreaGblForm(position, numVtx, v, basisByArea, gradBasisByArea);
            getBasisByAreaGblForm(position, numVtx, v, basisByArea);

            Vec2d wp_coord(0.0,0.0);
            double wp_grad = 0.0;
            for(int i=0; i< numVtx; i++){
                wp_coord = wp_coord + v[i]*basisByArea[i];
                wp_grad = wp_grad + gradBasisByArea[i].dot(v[i]);
            }
            assert(wp_coord[0] - MPsPosition(mp,0) < TEST_EPSILON);
            assert(wp_coord[1] - MPsPosition(mp,1) < TEST_EPSILON);
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
            for(int i=0; i<= numVtx; i++){
                wp_coord = wp_coord + v3d[i]*basisByArea3d[i];
                wp_coord2 = wp_coord + v3d[i]*basisByArea3d2[i];
            }
            /*printf("interpolation:(%.16e %.16e %.16e)\noriginal MP:(%.16e %.16e %.16e)\n",wp_coord[0],
                                              wp_coord[1],
                                              wp_coord[2],
                                              MPsPosition(mp,0),
                                              MPsPosition(mp,1),
                                              MPsPosition(mp,2));*/
            //assert(wp_coord[0] - MPsPosition(mp,0) < TEST_EPSILON);
            //assert(wp_coord[1] - MPsPosition(mp,1) < TEST_EPSILON);
            //assert(wp_coord[2] - MPsPosition(mp,2) < TEST_EPSILON);
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
            double basisByArea[maxVtxsPerElm] = {0.0};
            initArray(basisByArea,maxVtxsPerElm,0.0);
            Vec3d position(MPsPosition(mp,0),MPsPosition(mp,1),MPsPosition(mp,2));
            getBasisByAreaGblForm3d(position, numVtx, v, basisByArea);

            Vec3d wp_coord(0.0,0.0,0.0);
            for(int i=0; i< numVtx; i++){
                wp_coord = wp_coord + v[i]*basisByArea[i];
            }
            assert(wp_coord[0] - MPsPosition(mp,0) < TEST_EPSILON);
            assert(wp_coord[1] - MPsPosition(mp,1) < TEST_EPSILON);
	    assert(wp_coord[2] - MPsPosition(mp,2) < TEST_EPSILON);
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
