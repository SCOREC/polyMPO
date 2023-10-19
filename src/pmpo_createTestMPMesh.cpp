#include <Kokkos_Random.hpp>
#include <vector>
#include "pmpo_createTestMPMesh.hpp"

namespace polyMPO {

const int randSeed = 12345;

Mesh* createMesh(const mesh_type meshType, const geom_type geomType,
                 const double sphereRadius,
                 const int nVertices_size, const int nCells_size,
                 const std::vector<std::vector<double>> &v_array,
                 const std::vector<std::vector<int>> &elm2VtxConn_array,
                 const std::vector<int> &vtxCoords_array,
                 const std::vector<std::vector<int>> &elm2ElmConn_array){
    const int nCells = nCells_size;
    const int nVertices = nVertices_size;
    DoubleVec3dView vtxCoords("verticesCoordinates", nVertices);
    DoubleVec3dView::HostMirror h_vtxCoords = Kokkos::create_mirror_view(vtxCoords);
    for(int i=0; i<nVertices_size; i++){
        h_vtxCoords(i,0) = v_array[i][0]; 
        h_vtxCoords(i,1) = v_array[i][1];
        h_vtxCoords(i,2) = v_array[i][2]; 
    }
    Kokkos::deep_copy(vtxCoords, h_vtxCoords);

    IntVtx2ElmView elm2VtxConn("elementToVerticesConnection",nCells);
    IntVtx2ElmView::HostMirror h_elm2VtxConn = Kokkos::create_mirror_view(elm2VtxConn);
    for(int i=0; i<nCells_size; i++){
        h_elm2VtxConn(i,0) = vtxCoords_array[i];
        for(int j=0; j<h_elm2VtxConn(i,0); j++){
            h_elm2VtxConn(i,j+1) = elm2VtxConn_array[i][j];
        }
    }
    Kokkos::deep_copy(elm2VtxConn, h_elm2VtxConn);

    IntElm2ElmView elm2ElmConn("elementToElementsConnection",nCells);
    IntElm2ElmView::HostMirror h_elm2ElmConn = Kokkos::create_mirror_view(elm2ElmConn);
    for(int i=0; i<nCells_size; i++){
        h_elm2ElmConn(i,0) = vtxCoords_array[i];
        for(int j=0; j<h_elm2ElmConn(i,0); j++){
            h_elm2ElmConn(i,j+1) = elm2ElmConn_array[i][j];
        }
    }
    Kokkos::deep_copy(elm2ElmConn, h_elm2ElmConn);

    return new Mesh(meshType,
                    geomType,
                    sphereRadius,
                    nVertices,
                    nCells,
                    vtxCoords,
                    elm2VtxConn,
                    elm2ElmConn);
}

Mesh* initTestMesh(const int testMeshOption, const int replicateFactor){
    //add a switch statement testMeshOption
   if(testMeshOption==1 || testMeshOption==2){//this is a 2d mesh with 10 tri-oct elements
        const mesh_type meshType = mesh_general_polygonal;
        const geom_type geomType = geom_planar_surf;
        const double sphereRadius = 0.0;
        const int nVertices_size = 19;
        const int nCells_size = 10;
        std::vector<std::vector<double>> v_array_coords = //[nVertices_size][vec3d_nEntries]
                    {{0.00,0.00,1.1},{0.47,0.00,1.1},{1.00,0.00,1.1},{0.60,0.25,1.1},{0.60,0.40,1.1}, //5
                    {0.00,0.50,1.1},{0.31,0.60,1.1},{0.40,0.60,1.1},{0.45,0.55,1.1},{0.70,0.49,1.1}, //10
                    {0.80,0.45,1.1},{0.90,0.47,1.1},{1.00,0.55,1.1},{0.60,0.60,1.1},{0.37,0.80,1.1}, //15
                    {0.00,1.00,1.1},{0.37,1.00,1.1},{0.48,1.00,1.1},{1.00,1.00,1.1}};                //19
        
        if(testMeshOption==2){ 
            // rotation around x axis by 30 degrees
            double rad  = M_PI / 6.00;
            std::array<std::array<double, 3>, 3> R1 = {{{1.00, 0.00, 0.00}, 
                                                       {0.00, cos(rad), -sin(rad)},
                                                       {0.00, sin(rad), cos(rad) }}};
            // rotation around z axis by 45 degrees
            double rad2 = M_PI / 4.00;
            std::array<std::array<double, 3>, 3> R2 = {{{cos(rad2), -sin(rad2), 0.00}, 
                                                       {sin(rad2), cos(rad2), 0.00},
                                                       {0.00, 0.00, 1.00}}};

            for (int i = 0; i < (int)  v_array_coords.size(); i++) {
                std::array<double, 3> tmp_vec = {{0,0,0}};
                tmp_vec[0] = v_array_coords[i][0] * R1[0][0]
                           + v_array_coords[i][1] * R1[0][1]
                           + v_array_coords[i][2] * R1[0][2];
                tmp_vec[1] = v_array_coords[i][0] * R1[1][0]
                           + v_array_coords[i][1] * R1[1][1]
                           + v_array_coords[i][2] * R1[1][2];
                tmp_vec[2] = v_array_coords[i][0] * R1[2][0]
                           + v_array_coords[i][1] * R1[2][1]
                           + v_array_coords[i][2] * R1[2][2];

                v_array_coords[i][0] = tmp_vec[0] * R2[0][0]
                                     + tmp_vec[1] * R2[0][1]
                                     + tmp_vec[2] * R2[0][2];
                v_array_coords[i][1] = tmp_vec[0] * R2[1][0]
                                     + tmp_vec[1] * R2[1][1]
                                     + tmp_vec[2] * R2[1][2];
                v_array_coords[i][2] = tmp_vec[0] * R2[2][0]
                                     + tmp_vec[1] * R2[2][1]
                                     + tmp_vec[2] * R2[2][2];
            }
        }
        const std::vector<std::vector<double>> v_array = v_array_coords;
        const std::vector<std::vector<int>> elm2VtxConn_array = //[nCells_size][maxVtxsPerElm]
            {{1,2,4,5,9,8,7,6},       {2,3,4,-1,-1,-1,-1,-1},
             {4,3,11,10,5,-1,-1,-1},  {3,12,11,-1,-1,-1,-1,-1},
             {3,13,12,-1,-1,-1,-1,-1},{5,10,14,9,-1,-1,-1,-1},
             {9,14,18,17,15,8,-1,-1}, {7,8,15,-1,-1,-1,-1,-1},
             {6,7,15,17,16,-1,-1,-1}, {14,10,11,12,13,19,18,-1}};
        const std::vector<int> vtxCoords_array = {8,3,5,3,3,4,6,3,5,7}; //[nCells_size]
        const std::vector<std::vector<int>> elm2ElmConn_array = //[nCells_size][maxVtxsPerElm]
            {{-2,1,2,5,6,7,8,-2},    {-2,2,0,-1,-1,-1,-1,-1},
             {1,3,9,5,0,-1,-1,-1},   {4,9,2,-1,-1,-1,-1,-1},
             {-2,9,3,-1,-1,-1,-1,-1},{2,9,6,0,-1,-1,-1,-1},
             {5,9,-2,8,7,0,-1,-1},   {0,6,8,-1,-1,-1,-1,-1},
             {0,7,6,-2,-2,-1,-1,-1}, {5,2,3,4,-2,-2,6,-1}}; 
        auto meshReturn = createMesh(meshType, geomType, sphereRadius,
                                     nVertices_size, nCells_size,
                                     v_array, elm2VtxConn_array,
                                     vtxCoords_array,
                                     elm2ElmConn_array);
        if(replicateFactor >1){
            auto meshReplicate = replicateMesh(meshReturn,
                                               replicateFactor);
            delete meshReturn;
            return meshReplicate;
        }else{
            return meshReturn; 
        }
    }else{
        fprintf(stderr,"TestMeshOption not avaiable! return an empty mesh!");
        return new Mesh();
    }
    
}

MaterialPoints* initTestMPs(Mesh* mesh, int testMPOption){
    std::vector<int> mpPerElement;
    int numElms = mesh->getNumElements();
    switch(testMPOption){
        case 1:
            for(int i=0; i<numElms; i++){
                mpPerElement.push_back(i%3+4);
            }
        break;
        default:
            fprintf(stderr,"TestMPOption not avaiable! return an empty one!");
    }
    DoubleVec3dView vtxCoords = mesh->getMeshField<polyMPO::MeshF_VtxCoords>();   
    IntVtx2ElmView elm2VtxConn = mesh->getElm2VtxConn();
    auto geomType = mesh->getGeomType();

    int numMPs = 0;
    IntView numMPsPerElement("numMaterialPointsPerElement",numElms);
    if(mpPerElement.empty()) {
      Kokkos::Random_XorShift64_Pool<> random_pool(randSeed);
      Kokkos::parallel_for("setNumMPPerElement",numElms, KOKKOS_LAMBDA(const int i){
          auto generator = random_pool.get_state();
          numMPsPerElement(i) = generator.urand(4,7);
          random_pool.free_state(generator);
      });
    } else {
      PMT_ALWAYS_ASSERT(static_cast<size_t>(mesh->getNumElements()) == mpPerElement.size());
      Kokkos::View<int*, Kokkos::HostSpace> mpPerElement_hv(mpPerElement.data(),mpPerElement.size());
      Kokkos::deep_copy(numMPsPerElement,mpPerElement_hv);
    }

    Kokkos::parallel_reduce("calcTotalMP",numElms,KOKKOS_LAMBDA(const int&i, int& sum){
        sum += numMPsPerElement(i);
    },numMPs);
    IntView MPToElement("MPToElement",numMPs);

    Kokkos::parallel_scan("setMPsToElement", numElms, KOKKOS_LAMBDA(int i, int& iMP, bool is_final){
        if(is_final){  
            for(int j=0; j<numMPsPerElement(i); j++){
                MPToElement(iMP+j) = i;
            }   
        }
        iMP += numMPsPerElement(i); 
    },numMPs);

    DoubleVec3dView positions("MPpositions",numMPs);
    if(geomType == geom_planar_surf){     
        Kokkos::parallel_for("intializeMPsPositionPlanar", numMPs, KOKKOS_LAMBDA(const int iMP){
            int ielm = MPToElement(iMP);
            int numVtx = elm2VtxConn(ielm,0);
            double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
            for(int i=1; i<= numVtx; i++){
                sum_x += vtxCoords(elm2VtxConn(ielm,i)-1,0);
                sum_y += vtxCoords(elm2VtxConn(ielm,i)-1,1);
                sum_z += vtxCoords(elm2VtxConn(ielm,i)-1,2);
            }
            positions(iMP,0) = sum_x/numVtx;
            positions(iMP,1) = sum_y/numVtx;
            positions(iMP,2) = sum_z/numVtx;
        });
    }else if(geomType == geom_spherical_surf){
        Kokkos::Random_XorShift64_Pool<> random_pool(randSeed);
        const double radius = mesh->getSphereRadius();
        Kokkos::parallel_for("intializeMPsPositionSpherical", numMPs, KOKKOS_LAMBDA(const int iMP){
            int ielm = MPToElement(iMP);
            int numVtx = elm2VtxConn(ielm,0);
            auto generator = random_pool.get_state();
            int r1Index = generator.urand(0,numVtx);
            // int r2Index = (r1Index+numVtx/2)%numVtx;
            int r2Index = (r1Index+1)%numVtx;
            ++r1Index;
            ++r2Index;
            random_pool.free_state(generator);
            Vec3d midPoint(0.0,0.0,0.0);
            midPoint[0] = (vtxCoords(elm2VtxConn(ielm,r1Index)-1,0)+
                           vtxCoords(elm2VtxConn(ielm,r2Index)-1,0))/2;
            midPoint[1] = (vtxCoords(elm2VtxConn(ielm,r1Index)-1,1)+
                           vtxCoords(elm2VtxConn(ielm,r2Index)-1,1))/2;
            midPoint[2] = (vtxCoords(elm2VtxConn(ielm,r1Index)-1,2)+
                           vtxCoords(elm2VtxConn(ielm,r2Index)-1,2))/2;
            midPoint = midPoint * (1/midPoint.magnitude()) * radius; 
            positions(iMP,0) = midPoint[0];
            positions(iMP,1) = midPoint[1];
            positions(iMP,2) = midPoint[2];
        });
    } else{
        fprintf(stderr,"The geom type is not correct!");
        exit(1);
    }
    return new MaterialPoints(numElms,numMPs,positions,numMPsPerElement,MPToElement);
}

MPMesh initTestMPMesh(Mesh* mesh, int setMPOption) {
  auto mps = initTestMPs(mesh, setMPOption);
  return MPMesh(mesh,mps);
}

Mesh* replicateMesh(Mesh* mesh, int replicateFactor){
    //read the mesh
    const int nVertices_size = mesh->getNumVertices();
    const int nCells_size = mesh->getNumElements();
    const mesh_type meshType = mesh->getMeshType(); 
    const geom_type geomType = mesh->getGeomType(); 
    const double sphereRadius = mesh->getSphereRadius(); 
    const int nCells = nCells_size*replicateFactor;
    const int nVertices = nVertices_size*replicateFactor;
     
    DoubleVec3dView v_array = mesh->getMeshField<polyMPO::MeshF_VtxCoords>(); 
    DoubleVec3dView vtxCoords("verticesCoordinates", nVertices);
    Kokkos::parallel_for("set vtxCoords", nVertices_size, KOKKOS_LAMBDA(const int i){
        for(int f=0; f<replicateFactor; f++){
            vtxCoords(i+f*nVertices_size,0) = v_array(i,0); 
            vtxCoords(i+f*nVertices_size,1) = v_array(i,1);
            vtxCoords(i+f*nVertices_size,2) = v_array(i,2); 
        }
    });

    IntVtx2ElmView elm2VtxConn_array = mesh->getElm2VtxConn();
    IntVtx2ElmView elm2VtxConn("elementToVerticesConnection",nCells);
    Kokkos::parallel_for("set elm2VtxConn", nCells_size, KOKKOS_LAMBDA(const int i){
        for(int f=0; f<replicateFactor; f++){
            elm2VtxConn(i+f*nCells_size,0) = elm2VtxConn_array(i,0);
            for(int j=0; j<elm2VtxConn(i,0); j++){
                elm2VtxConn(i+f*nCells_size,j+1) = elm2VtxConn_array(i,j+1) + f*nVertices_size;
            }
        }
    });

    IntElm2ElmView elm2ElmConn_array = mesh->getElm2ElmConn();
    IntElm2ElmView elm2ElmConn("elementToElementsConnection",nCells);
    Kokkos::parallel_for("set elm2ElmConn", nCells_size, KOKKOS_LAMBDA(const int i){
        for(int f=0; f<replicateFactor; f++){
            elm2ElmConn(i+f*nCells_size,0) = elm2ElmConn_array(i,0);
            for(int j=0; j<elm2ElmConn(i,0); j++){
                elm2ElmConn(i+f*nCells_size,j+1) = elm2ElmConn_array(i,j+1) + f*nVertices_size;
            }
        }
    });
    
    //create new mesh   
    return new Mesh(meshType,
                    geomType,
                    sphereRadius,
                    nVertices,
                    nCells,
                    vtxCoords,
                    elm2VtxConn,
                    elm2ElmConn);

}

}
