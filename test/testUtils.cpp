#include <Kokkos_Random.hpp>
#include "testUtils.hpp"

Mesh initTestMesh(int factor){
    const int nCells_size = 10;//10 tri-oct
    const int nVertices_size = 19;//arbi <20
    int nCells = nCells_size*factor;
    int nVertices = nVertices_size*factor;
    const double v_array[19][2] = {{0.00,0.00},{0.47,0.00},{1.00,0.00},{0.60,0.25},{0.60,0.40},{0.00,0.50},{0.31,0.60},{0.40,0.60},{0.45,0.55},{0.70,0.49},{0.80,0.45},{0.90,0.47},{1.00,0.55},{0.60,0.60},{0.37,0.80},{0.00,1.00},{0.37,1.00},{0.48,1.00},{1.00,1.00}};
    const int elm2VtxConn_array[10][8] = {{1,2,4,5,9,8,7,6},{2,3,4,-1,-1,-1,-1,-1},
        {4,3,11,10,5,-1,-1,-1},{3,12,11,-1,-1,-1,-1,-1},{3,13,12,-1,-1,-1,-1,-1},
        {5,10,14,9,-1,-1,-1,-1},{9,14,18,17,15,8,-1,-1},{7,8,15,-1,-1,-1,-1,-1},
        {6,7,15,17,16,-1,-1,-1},{14,10,11,12,13,19,18,-1}};
    const int vtxCoords_array[10] = {8,3,5,3,3,4,6,3,5,7};
    const int vtx2ElmConn_array[19][6] = {{1,0,-1,-1,-1,-1},{2,0,1,-1,-1,-1},{4,1,2,3,4,-1},
        {3,0,1,2,-1,-1},{3,0,2,5,-1,-1},{2,0,8,-1,-1,-1},{3,0,7,8,-1,-1},
        {3,0,6,7,-1,-1},{3,0,5,6,-1,-1},{3,2,5,9,-1,-1},{3,2,3,9,-1,-1},
        {3,3,4,9,-1,-1},{2,4,9,-1,-1,-1},{3,5,6,9,-1,-1},{3,6,7,8,-1,-1},
        {1,8,-1,-1,-1,-1},{2,6,8,-1,-1,-1},{2,6,9,-1,-1,-1},{1,9,-1,-1,-1,-1}};
    //-1 is null, -2 is not connected to any elm
    const int elm2ElmConn_array[10][8] = {{-2,1,2,5,6,7,8,-2},{-2,2,0,-1,-1,-1,-1,-1},
        {1,3,9,5,0,-1,-1,-1},{4,9,2,-1,-1,-1,-1,-1},{-2,9,3,-1,-1,-1,-1,-1},
        {2,9,6,0,-1,-1,-1,-1},{5,9,-2,8,7,0,-1,-1},{0,6,8,-1,-1,-1,-1,-1},
        {0,7,6,-2,-2,-1,-1,-1},{5,2,3,4,-2,-2,6,-1}}; 
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
    
    IntElm2ElmView elm2ElmConn("elementToElementsConnection",nCells);
    IntElm2ElmView::HostMirror h_elm2ElmConn = Kokkos::create_mirror_view(elm2ElmConn);
    for(int f=0; f<factor; f++){
        for(int i=0; i<nCells_size; i++){
            h_elm2ElmConn(i+f*nCells_size,0) = vtxCoords_array[i];
            for(int j=0; j<h_elm2ElmConn(i,0); j++){
                h_elm2ElmConn(i+f*nCells_size,j+1) = elm2ElmConn_array[i][j] + f*nVertices_size;
            }
        }
    }
    Kokkos::deep_copy(elm2ElmConn, h_elm2ElmConn);

    return Mesh(nVertices,
                nCells,
                vtxCoords,
                elm2VtxConn,
                vtx2ElmConn,
                elm2ElmConn);
}


MPM initTestMPM(Mesh& mesh){
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    IntVtx2ElmView elm2VtxConn = mesh.getElm2VtxConn();

    int numMPs = 0;
    
    IntView numMPsPerElement("numMaterialPointsPerElement",numElms);
    Kokkos::Random_XorShift64_Pool<> random_pool(randSeed);
    Kokkos::parallel_for("setNumMPPerElement",numElms, KOKKOS_LAMBDA(const int i){
        auto generator = random_pool.get_state();
        numMPsPerElement(i) = generator.urand(4,7);   
        random_pool.free_state(generator);
    });
    Kokkos::fence();

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
    });
    
    auto p = MaterialPoints(numMPs,positions,isActive);
    return MPM(mesh,p,elm2MPs,MPs2Elm);
}


void interpolateWachspress(MPM& mpm){
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();   
    auto MPs2Elm = mpm.getMPs2Elm();
   
    auto vtxCoords = mesh.getVtxCoords();
    auto elm2VtxConn = mesh.getElm2VtxConn();

    auto numMPs = MPs.getCount();
    auto MPsPosition = MPs.getPositions();
    auto isActive = MPs.isActive();

    Kokkos::parallel_for("getBasisAndGradByAreaGblForm",numMPs,KOKKOS_LAMBDA(const int iMP){
        if (isActive(iMP)){
            int iElm = MPs2Elm(iMP);
            Vector2 v[maxVtxsPerElm+1] = {vtxCoords(elm2VtxConn(iElm,1))};
            initArray(v,maxVtxsPerElm+1,vtxCoords(elm2VtxConn(iElm,1)));
            int numVtx = elm2VtxConn(iElm,0);
            for(int i = 1; i<=numVtx; i++){
                v[i-1] = vtxCoords(elm2VtxConn(iElm,i)-1);
            }
            v[numVtx] = vtxCoords(elm2VtxConn(iElm,1)-1);
            double basisByArea[maxVtxsPerElm] = {0.0};
            initArray(basisByArea,maxVtxsPerElm,0.0);
            Vector2 gradBasisByArea[maxVtxsPerElm];
            getBasisAndGradByAreaGblForm(MPsPosition(iMP), numVtx, v, basisByArea, gradBasisByArea);
        
            Vector2 wp_coord(0.0,0.0);
            double wp_grad = 0.0;
            for(int i=0; i<= numVtx; i++){
                wp_coord = wp_coord + v[i]*basisByArea[i];
                wp_grad = wp_grad + gradBasisByArea[i].dot(v[i]);
            }
        }        
    });
}

Vector2View initT2LDelta(int size, const double range, const int randomSeed){
    Vector2View retVal("T2LDeltaX",size);

    Kokkos::Random_XorShift64_Pool<> random_pool(randomSeed);
    Kokkos::parallel_for("setT2LDeltaX", size, KOKKOS_LAMBDA(const int i){
        auto generator = random_pool.get_state();
        retVal(i) = Vector2(generator.drand(-range,range),generator.drand(-range,range));
        random_pool.free_state(generator);
    });
    Kokkos::fence();

    return retVal;
}


MPM initMPMWithRandomMPs(Mesh& mesh, int factor, const int randomSeed){
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    IntVtx2ElmView elm2VtxConn = mesh.getElm2VtxConn(); 

    int numMPs = numElms*factor;
    Vector2View positions("MPpositions", numMPs);
    BoolView isActive("MPstatus",numMPs);
    IntView elm2MPs("elementToMPs",numElms*(maxMPsPerElm+1));
    IntView MPs2Elm("MPToElementIDs",numMPs);

    IntView MPToElement("MPToElement",numMPs);
    Kokkos::Random_XorShift64_Pool<> random_pool(randomSeed);

     Kokkos::parallel_scan("setMPsToElement", numElms, KOKKOS_LAMBDA(int i, int& iMP, bool is_final){
        if(is_final){  
            elm2MPs(i*(maxMPsPerElm+1)) = factor;
            for(int j=0; j<factor; j++){
                MPToElement(iMP+j) = i;
                elm2MPs(i*(maxMPsPerElm+1)+j+1) = iMP+j;
            }   
        }
        iMP += factor; 
    },numMPs);

    Kokkos::parallel_for("setPositions", numMPs, KOKKOS_LAMBDA(const int iMP){
        auto generator = random_pool.get_state();
        int ielm = MPToElement(iMP);
        int numVtx = elm2VtxConn(ielm,0);
        double sum_x = 0.0, sum_y = 0.0;
        for(int i=1; i<= numVtx; i++){
            sum_x += vtxCoords(elm2VtxConn(ielm,i)-1)[0];
            sum_y += vtxCoords(elm2VtxConn(ielm,i)-1)[1];
        }
        Vector2 XYc = Vector2(sum_x/numVtx, sum_y/numVtx);
        int triID = generator.urand(0,numVtx);
        double rws[2] = {generator.drand(0.0,1.0), generator.drand(0.0,1.0)};
        random_pool.free_state(generator);
        double weights[3];
        if (rws[0]> rws[1]){
            weights[0] = rws[1];
            weights[1] = rws[0]-rws[1];
            weights[2] = 1-rws[0];
        }else{
            weights[0] = rws[0];
            weights[1] = rws[1]-rws[0];
            weights[2] = 1-rws[1];
        }
        auto v1 = vtxCoords(elm2VtxConn(ielm,triID+1)-1);
        auto v2 = vtxCoords(elm2VtxConn(ielm,(triID+1)%numVtx+1)-1);
        positions(iMP) = XYc*weights[0]+v1*weights[1]+v2*weights[2];
        MPs2Elm(iMP) = ielm;
        isActive(iMP) = true;
    });    

    auto p = MaterialPoints(numMPs,positions,isActive);
    return MPM(mesh,p,elm2MPs,MPs2Elm);
}

Vector2View initT2LDeltaRankineVortex(MPM mpm, Vector2 center, const int numEdge, const double dx, const double Gamma ){
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();
    int numMPs = MPs.getCount();
    auto MPsPosition = MPs.getPositions();    

    Vector2View returnDx("T2LDeltaXY",numMPs);
    const double a = numEdge*dx;
    const double coeff = Gamma/(2*MPMTEST_PI);
    //T = (2*MPMTEST_PI*a)*(2*MPMTEST_PI*a)/Gamma;
    const double dt = 4*dx*2*MPMTEST_PI*a/Gamma; // a= numEdge*dx
    const double coeffLess = coeff*dt/(a*a);
    const double coeffGret = coeff*dt;
    Kokkos::parallel_for("setNumMPPerElement", numMPs, KOKKOS_LAMBDA(const int iMP){
        auto MP = MPsPosition(iMP);
        auto centerVector = MP-center;
        auto radius = centerVector.magnitude();
        double vTheta;
        if(radius<= a){
            vTheta = coeffLess*radius;
        }else{
            vTheta = coeffGret/radius;
        }
        //(-y,+x) to get tangential
        returnDx(iMP) = Vector2(-centerVector[1],centerVector[0])*(vTheta/radius);
    });
    Kokkos::fence();

    return returnDx;
}

double calcAvgLengthOfEdge(Mesh mesh, int printOut){
    auto vtxCoords = mesh.getVtxCoords();
    auto elm2VtxConn = mesh.getElm2VtxConn();
    auto elm2ElmConn = mesh.getElm2ElmConn();
    auto numElm = mesh.getNumElements();

    double sum = 0.0;
    double SqrSum = 0.0;
    int count = 0;
    int countTwice = 0;
    Kokkos::parallel_reduce("sumOfEdges*2", numElm,KOKKOS_LAMBDA(const int& iElm, double& lsum, int& lcount, double& lSqrSum, int& lcountTwice){
        int numVtx = elm2VtxConn(iElm,0);
        Vector2 v[maxVtxsPerElm+1];
        for(int i = 1; i<=numVtx; i++){
                v[i-1] = vtxCoords(elm2VtxConn(iElm,i)-1);
        }
        v[numVtx] = vtxCoords(elm2VtxConn(iElm,1)-1);
        for (int i = 0; i < numVtx; i++){
            lcount++;
            double l_e = (v[i+1] - v[i]).magnitude();
            if( elm2ElmConn(iElm,i+1)>= 0){//connected edge counted twice
                lsum += l_e/2;
                lSqrSum += l_e*l_e/2;
                lcountTwice++;
            }else{
                lsum += l_e;
                lSqrSum += l_e*l_e;
            }
        }
    },sum,count,SqrSum,countTwice);
    count -= countTwice/2;
    if(printOut)
        printf("%d: l_eSum= %f, avg_l_e= %f, avg_l_e*l_e= %f\n",count, sum, sum/count, SqrSum/count);
    return sum/count;
}

Vector2View initT2LTest1(const int size, const double range, double percent1, double percent2, double percent3, double percent4, const int randomSeed){
    PMT_ALWAYS_ASSERT(percent1+percent2+percent3+percent4 <= 1+MPMTEST_EPSILON);
    percent2 += percent1;
    percent3 += percent2;
    percent4 += percent3;
    Vector2View returnDx("T2LDeltaX",size);
    return returnDx;
}

Vector2View initT2LTest2(const int size, const double range, double percent1, double percent2, double percent3, double percent4, const int randomSeed){ 
    PMT_ALWAYS_ASSERT(percent1+percent2+percent3+percent4 <= 1+MPMTEST_EPSILON);
    percent2 += percent1;
    percent3 += percent2;
    percent4 += percent3;
    Vector2View returnDx("T2LDeltaX",size);    
    Kokkos::Random_XorShift64_Pool<> random_pool(randomSeed);
    Kokkos::parallel_for("setNumMPPerElement", size, KOKKOS_LAMBDA(const int i){
        auto generator = random_pool.get_state();
        double direction = generator.drand(0,2*MPMTEST_PI);
        double iRange = generator.drand(1);
        double length = 0.0;
        if(iRange < percent1){
            length = generator.drand(0.0,range*0.5);
        }else if(iRange < percent2){
            length = generator.drand(range,range*1.5);
        }else if(iRange < percent3){
            length = generator.drand(range*2,range*2.5);
        }else{
            length = generator.drand(range*3,range*3.5);
        }
        returnDx(i) = Vector2(length*std::cos(direction), length*std::sin(direction));
        random_pool.free_state(generator);
    });
    Kokkos::fence();
    
    return returnDx;
}
