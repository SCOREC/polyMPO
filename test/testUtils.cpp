#include <Kokkos_Random.hpp>
#include "testUtils.hpp"

const int randSeed = 12345;

Mesh initTestMesh(int factor){
    const int nCells_size = 10;//10 tri-oct
    const int nVertices_size = 19;//arbi <20
    int nCells = nCells_size*factor;
    int nVertices = nVertices_size*factor;
    const double v_array[19][2] = {{0.00,0.00},{0.47,0.00},{1.00,0.00},{0.60,0.25},{0.60,0.40},{0.00,0.50},{0.31,0.60},{0.40,0.60},{0.45,0.55},{0.70,0.49},{0.80,0.45},{0.90,0.47},{1.00,0.55},{0.60,0.60},{0.37,0.80},{0.00,1.00},{0.37,1.00},{0.48,1.00},{1.00,1.00}};
    const int elm2VtxConn_array[10][8] = {{1,2,4,5,9,8,7,6},{2,3,4,-1,-1,-1,-1,-1},{4,3,11,10,5,-1,-1,-1},
    {3,12,11,-1,-1,-1,-1,-1},{3,13,12,-1,-1,-1,-1,-1},
    {5,10,14,9,-1,-1,-1,-1},{9,14,18,17,15,8,-1,-1},
    {7,8,15,-1,-1,-1,-1,-1},{6,7,15,17,16,-1,-1,-1},
    {14,10,11,12,13,19,18,-1}};
    const int vtxCoords_array[10] = {8,3,5,3,3,4,6,3,5,7};
    const int vtx2ElmConn_array[19][6] = {{1,0,-1,-1,-1,-1},{2,0,1,-1,-1,-1},{4,1,2,3,4,-1},{3,0,1,2,-1,-1},{3,0,2,5,-1,-1},{2,0,8,-1,-1,-1},{3,0,7,8,-1,-1},{3,0,6,7,-1,-1},{3,0,5,6,-1,-1},{3,2,5,9,-1,-1},{3,2,3,9,-1,-1},{3,3,4,9,-1,-1},{2,4,9,-1,-1,-1},{3,5,6,9,-1,-1},{3,6,7,8,-1,-1},{1,8,-1,-1,-1,-1},{2,6,8,-1,-1,-1},{2,6,9,-1,-1,-1},{1,9,-1,-1,-1,-1}};
    
    Vector2View vtxCoords("verticesCoordinates", nVertices);
    IntElm2VtxView vtx2ElmConn("vertexToElementsConnection",nVertices);
    
    Vector2View::HostMirror h_vtxCoords = Kokkos::create_mirror_view(vtxCoords);
    IntElm2VtxView::HostMirror h_vtx2ElmConn = Kokkos::create_mirror_view(vtx2ElmConn);
    for(int f=0; f<factor; f++){
        for(int i=0; i<nVertices_size; i++){
            h_vtxCoords(i+f*nVertices_size) = Vector2(v_array[i][0],v_array[i][1]); 
            h_vtx2ElmConn(i+f*nVertices_size,0) = vtx2ElmConn_array[i][0];
            for(int j=1; j<=vtx2ElmConn_array[i][0]; j++){
                h_vtx2ElmConn(i+f*nVertices_size,j) = vtx2ElmConn_array[i][j]+f*nCells_size; 
            }
        }
    }
    Kokkos::deep_copy(vtxCoords, h_vtxCoords);
    Kokkos::deep_copy(vtx2ElmConn,h_vtx2ElmConn);

    IntVtx2ElmView elm2VtxConn("elementToVerticesConnection",nCells);
    
    IntVtx2ElmView::HostMirror h_elm2VtxConn = Kokkos::create_mirror_view(elm2VtxConn);
    for(int f=0; f<factor; f++){
        for(int i=0; i<nCells_size; i++){
            h_elm2VtxConn(i+f*nCells_size,0) = vtxCoords_array[i];
            for(int j=0; j<h_elm2VtxConn(i,0); j++){
                h_elm2VtxConn(i+f*nCells_size,j+1) = elm2VtxConn_array[i][j] + f*nVertices_size;
            }
        }
    }
    Kokkos::deep_copy(elm2VtxConn, h_elm2VtxConn);
    
    return Mesh(nVertices,
                nCells,
                vtxCoords,
                elm2VtxConn,
                vtx2ElmConn);
}


MPMesh initTestMPMesh(Mesh& mesh){
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    IntVtx2ElmView elm2VtxConn = mesh.getElm2VtxConn();

    int numMPs = 0;
    
    IntView numMPsPerElement("numMaterialPointsPerElement",numElms);
    Kokkos::Random_XorShift64_Pool<> random_pool(randSeed);
    Kokkos::parallel_for("setNumMPPerElement",numElms, KOKKOS_LAMBDA(const int i){
        auto generator = random_pool.get_state();
        numMPsPerElement(i) = generator.urand(4,7); //rand between 4 and 7 - TODO: make input arg
        random_pool.free_state(generator);
    });

    Kokkos::parallel_reduce("calcTotalMP",numElms,KOKKOS_LAMBDA(const int&i, int& sum){
        sum += numMPsPerElement(i);
    },numMPs);
    IntView MPToElement("MPToElement",numMPs);  

    IntView elm2MPs("elementToMPs",numElms*(maxMPsPerElm+1));
    Kokkos::parallel_scan("setMPsToElement", numElms, KOKKOS_LAMBDA(int i, int& iMP, bool is_final){
        if(is_final){  
            elm2MPs(i*(maxMPsPerElm+1)) = numMPsPerElement(i);
            for(int j=0; j<numMPsPerElement(i); j++){
                MPToElement(iMP+j) = i;
                elm2MPs(i*(maxMPsPerElm+1)+j+1) = iMP+j;
            }   
        }
        iMP += numMPsPerElement(i); 
    },numMPs);

    Vector2View positions("MPpositions",numMPs);
    IntView MPs2Elm("MPToElementIDs",numMPs);
    BoolView isActive("MPstatus",numMPs);
     
    Kokkos::parallel_for("intializeMPsPosition", numMPs, KOKKOS_LAMBDA(const int iMP){
        int ielm = MPToElement(iMP);
        int numVtx = elm2VtxConn(ielm,0);
        double sum_x = 0.0, sum_y = 0.0;
        for(int i=1; i<= numVtx; i++){
            sum_x += vtxCoords(elm2VtxConn(ielm,i)-1)[0];
            sum_y += vtxCoords(elm2VtxConn(ielm,i)-1)[1];
        }
        positions(iMP) = Vector2(sum_x/numVtx, sum_y/numVtx);
        MPs2Elm(iMP) = ielm;
        isActive(iMP) = true;
        //printf("%d: (%f,%f)\n",iMP,positions(iMP)[0],positions(iMP)[1]);
    });
    
    auto p = new MaterialPoints(numElms,numMPs,positions,isActive);

    return MPMesh(mesh,p,elm2MPs,MPs2Elm);
}


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
            getBasisAndGradByAreaGblForm(MPsPosition(mp), numVtx, v, basisByArea, gradBasisByArea);
        
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
