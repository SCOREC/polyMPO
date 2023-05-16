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
    auto elm2VtxConn = mesh.getElm2VtxConn();

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
        int iElm = MPToElement(iMP);
        positions(iMP) = calcElmCenter(iElm,elm2VtxConn,vtxCoords);
        MPs2Elm(iMP) = iElm;
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
    PMT_ALWAYS_ASSERT(range > 0);
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
    auto elm2VtxConn = mesh.getElm2VtxConn(); 

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
        int iElm = MPToElement(iMP);
        int numVtx = elm2VtxConn(iElm,0);
        Vector2 XYc = calcElmCenter(iElm,elm2VtxConn,vtxCoords);
//        if(ielm == 8 || ielm == 2420 || ielm == 4880 || ielm == 5286)
//            printf("%f %f 0.0\n", XYc[0], XYc[1]);
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
        auto v1 = vtxCoords(elm2VtxConn(iElm,triID+1)-1);
        auto v2 = vtxCoords(elm2VtxConn(iElm,(triID+1)%numVtx+1)-1);
        positions(iMP) = XYc*weights[0]+v1*weights[1]+v2*weights[2];
        MPs2Elm(iMP) = iElm;
        isActive(iMP) = true;
    });    

    auto p = MaterialPoints(numMPs,positions,isActive);
    return MPM(mesh,p,elm2MPs,MPs2Elm);
}

Vector2View initT2LDeltaRankineVortex(MPM mpm, Vector2 center, const int numEdge, const double dx, const double Gamma ){
    PMT_ALWAYS_ASSERT(numEdge >0);
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

Vector2View initT2LTest1(MPM mpm, double ratio1, double ratio2, double ratio3, double ratio4, const int randomSeed){
    PMT_ALWAYS_ASSERT(ratio1+ratio2+ratio3+ratio4 <= 1+MPMTEST_EPSILON);
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    int numMPs = MPs.getCount();
    auto MPsPosition = MPs.getPositions();
    auto MPs2Elm = mpm.getMPs2Elm();
    auto elm2ElmConn = mesh.getElm2ElmConn();
    auto isActive = MPs.isActive();
    auto elm2VtxConn = mesh.getElm2VtxConn();
    Vector2View returnDx("T2LDeltaX",numMPs);
    Kokkos::Random_XorShift64_Pool<> random_pool(randomSeed);
    //Ind = 1 bound; Ind = 2 connect to bound
    IntView MPBoundInd("MPatBoundaryIndicator",numMPs);
    IntView ElmBoundInd("MPatBoundaryIndicator",elm2ElmConn.size());
    IntView countBounds("countBounds",2);
    IntView::HostMirror h_countBounds = Kokkos::create_mirror_view(countBounds);
    Kokkos::parallel_for("setBoundary1", numMPs, KOKKOS_LAMBDA(const int iMP){
    if(isActive(iMP)){
        int iElm = MPs2Elm(iMP);
        int numEdge = elm2ElmConn(iElm, 0);
        for(int i=0; i<numEdge; i++){
            if(elm2ElmConn(iElm,i+1)<0){
                MPBoundInd(iMP) = 1;
                Kokkos::atomic_increment(&countBounds(0));
                ElmBoundInd(iElm) = 1;
                break;
            }
        }
    }});
    Kokkos::fence();
    Kokkos::parallel_for("setBoundary2", numMPs, KOKKOS_LAMBDA(const int iMP){
    if(isActive(iMP) && MPBoundInd(iMP)!= 1){
        int iElm = MPs2Elm(iMP);
        int numEdge = elm2ElmConn(iElm, 0);
        for(int i=0; i<numEdge; i++){
            if(ElmBoundInd(elm2ElmConn(iElm,i+1)) == 1){
                MPBoundInd(iMP) = 2;
                Kokkos::atomic_increment(&countBounds(1));
                ElmBoundInd(iElm) = 2;
                break;
            }
        }
    }});
    Kokkos::fence();
    deep_copy(h_countBounds, countBounds);
    int leftFreeMPs = numMPs - h_countBounds(0) - h_countBounds(1);
    double percent1 = (numMPs*ratio1-h_countBounds(0))/leftFreeMPs;
    double percent2 = (numMPs*ratio2-h_countBounds(1))/leftFreeMPs;
    double percent3 = (numMPs*ratio3)/leftFreeMPs;
    double percent4 = (numMPs*ratio4)/leftFreeMPs;
    double percentForPartOf2 = 1.0;
    if(percent2<0){
        if(ratio2 == 0.0){
            percent1 = (numMPs*ratio1-h_countBounds(0)-h_countBounds(1))/leftFreeMPs;
            percent2 = 0;
            percentForPartOf2 = 0.0;
        }else{
            percent1 = (numMPs*ratio1-h_countBounds(0)-h_countBounds(1)+numMPs*ratio2)/leftFreeMPs;
            percent2 = 0;
            percentForPartOf2 = numMPs*ratio2/h_countBounds(1);
        }
    }
    percent2 += percent1;
    percent3 += percent2;
    percent4 += percent3;
    Kokkos::parallel_for("setNumMPPerElement", numMPs, KOKKOS_LAMBDA(const int iMP){
    if(isActive(iMP)){
        auto generator = random_pool.get_state();
        int numAcross = 0;
        if(MPBoundInd(iMP) == 1){
            numAcross = 0;
        }else if(MPBoundInd(iMP) == 2){
            double range2 = generator.drand(1);
            if(range2 < percentForPartOf2){
                numAcross = 1;
            }else{
                numAcross = 0;
            }
        }else{
            double iRange = generator.drand(percent4);
            if(iRange < percent1){//maybe add a goto to improve a bit of performance
                numAcross = 0;
            }else if(iRange < percent2){
                numAcross = 1;
            }else if(iRange < percent3){
                numAcross = 2;
            }else if(iRange < percent4){
                numAcross = 3;
            }
        }
        int initElm = MPs2Elm(iMP);
        int iElm = initElm;
        int numEdge = elm2ElmConn(iElm, 0);
        int initEdge,nextElm;
        do{
            initEdge = generator.urand(0,numEdge)+1;
            nextElm = elm2ElmConn(iElm,initEdge); 
        }while(nextElm<0);
        random_pool.free_state(generator);
        Vector2 targetPosition = calcElmCenter(initElm,elm2VtxConn,vtxCoords);
        MPsPosition(iMP) = targetPosition;
        Vector2 MPPosition = MPsPosition(iMP);
        Vector2 v0 = vtxCoords(elm2VtxConn(initElm,initEdge)-1);
        Vector2 v1 = vtxCoords(elm2VtxConn(initElm,initEdge%numEdge+1)-1);
        Vector2 edgeCenter = (v1 + v0)*0.5;
        Vector2 direction = edgeCenter - targetPosition;
        direction = direction*(1/direction.magnitude());
        for(int iAcross = 0; iAcross< numAcross; iAcross++){
            if(nextElm <0){
                break;
            }
            int numVtx = elm2VtxConn(nextElm,0);
            int v[maxVtxsPerElm];
            for(int i=0; i< numVtx; i++)
                v[i] = elm2VtxConn(nextElm,i+1)-1;
            Vector2 e[maxVtxsPerElm];
            double pdx[maxVtxsPerElm];
            for(int i=0; i< numVtx; i++){
                Vector2 v_i = vtxCoords(v[i]);
                Vector2 v_ip1 = vtxCoords(v[(i+1)%numVtx]);
                e[i] = v_ip1 - v_i;
                pdx[i] = (v_i - targetPosition).cross(direction);
            }
            
            // update the nextElm and targetPosition  
            int tempForNextElm = initEdge;
            for(int i=0; i<numVtx; i++){
                int ip1 = (i+1)%numVtx;
                if(pdx[i]*pdx[ip1] <0){
                    if(elm2ElmConn(nextElm,i+1) == iElm){
                        Vector2 P2Vi = vtxCoords(v[i])-targetPosition; 
                        Vector2 edgeIntersect = targetPosition + direction*direction.dot(P2Vi);
                        targetPosition = edgeIntersect + direction*direction.dot(calcElmCenter(nextElm,elm2VtxConn,vtxCoords)-edgeIntersect);
                    }
                    else{
                        // i+1 will be the next edge to across
                        tempForNextElm = i+1;
                    }
                }
            }
            iElm = nextElm; 
            nextElm = elm2ElmConn(nextElm,tempForNextElm);
        }
        returnDx(iMP) = targetPosition - MPPosition;
    }});
    return returnDx;
}

Vector2View initT2LTest2(MPM mpm, const int MPAcross, const int randomSeed){ 
    PMT_ALWAYS_ASSERT(MPAcross >= 0);
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    int numMPs = MPs.getCount();
    auto MPsPosition = MPs.getPositions();
    auto MPs2Elm = mpm.getMPs2Elm();
    auto elm2ElmConn = mesh.getElm2ElmConn();
    auto isActive = MPs.isActive();
    IntVtx2ElmView elm2VtxConn = mesh.getElm2VtxConn();
    Vector2View returnDx("T2LDeltaX",numMPs);
    Kokkos::Random_XorShift64_Pool<> random_pool(randomSeed);
    Kokkos::parallel_for("setNumMPPerElement", numMPs, KOKKOS_LAMBDA(const int iMP){
    if(isActive(iMP)){
        auto generator = random_pool.get_state();
        Vector2 MPPosition = MPsPosition(iMP);
        int numAcross = 0;
        if(iMP == arbitraryInt){
            numAcross = MPAcross;
        }
        int initElm = MPs2Elm(iMP);
        int iElm = initElm;
        int numElm = elm2ElmConn(iElm, 0);
        int initDirection = generator.urand(0,numElm)+1;//[0:numElm)+1
        int nextEdge = initDirection;//[0:numElm)+1
        int nextElm = elm2ElmConn(iElm,nextEdge);
        for(int iAcross = 0; iAcross< numAcross; iAcross++){
            if(nextElm < 0){//reInit to a new direction
                iElm = initElm;
                numElm = elm2ElmConn(iElm, 0);
                if(initDirection == numElm){
                    initDirection = 1;
                }else{
                    initDirection += 1;
                }
                nextEdge = initDirection;
                nextElm = elm2ElmConn(iElm,nextEdge);
                iAcross = -1;
                continue;
            }
            int oldElm = iElm;
            iElm = nextElm;
            numElm = elm2ElmConn(iElm, 0);
            //update the nextEdge:
            for(int i=1; i<=numElm; i++){
                if(elm2ElmConn(iElm, i) == oldElm){
                    nextEdge = i;
                    break;
                }
            }
            nextEdge = (nextEdge + numElm/2)%numElm;
            if(nextEdge == 0)
                nextEdge = numElm;
            nextElm = elm2ElmConn(iElm,nextEdge);
        }
        int numVtx = elm2VtxConn(iElm,0);
        Vector2 XYc = calcElmCenter(iElm, elm2VtxConn, vtxCoords);
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
        auto v1 = vtxCoords(elm2VtxConn(iElm,triID+1)-1);
        auto v2 = vtxCoords(elm2VtxConn(iElm,(triID+1)%numVtx+1)-1);
        Vector2 targetPosition = XYc*weights[0]+v1*weights[1]+v2*weights[2];
        returnDx(iMP) = targetPosition - MPPosition;
    }}); 
    Kokkos::fence();
    
    return returnDx;
}

Vector2View initT2LTest3(MPM mpm, const int MPAcross, const double percent, const int randomSeed){ 
    PMT_ALWAYS_ASSERT(MPAcross >= 0);
    PMT_ALWAYS_ASSERT(percent >= 0 && percent <1.0+MPMTEST_EPSILON);
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    int numMPs = MPs.getCount();
    auto MPsPosition = MPs.getPositions();
    auto MPs2Elm = mpm.getMPs2Elm();
    auto elm2ElmConn = mesh.getElm2ElmConn();
    auto isActive = MPs.isActive();
    IntVtx2ElmView elm2VtxConn = mesh.getElm2VtxConn();
    Vector2View returnDx("T2LDeltaX",numMPs);
    Kokkos::Random_XorShift64_Pool<> random_pool(randomSeed);
    Kokkos::parallel_for("setNumMPPerElement", numMPs, KOKKOS_LAMBDA(const int iMP){
    if(isActive(iMP)){
        Vector2 MPPosition = MPsPosition(iMP);
        auto generator = random_pool.get_state();
        int numAcross = 0;
        if(generator.drand(1)<percent){
            numAcross = MPAcross;
        }
        random_pool.free_state(generator);
        int initElm = MPs2Elm(iMP);
        int iElm = initElm;
        int numElm = elm2ElmConn(iElm, 0);
        int initDirection = 1;//[0:numElm)+1
        int nextEdge = initDirection;//[0:numElm)+1
        int nextElm = elm2ElmConn(iElm,nextEdge);
        for(int iAcross = 0; iAcross< numAcross; iAcross++){
            if(nextElm < 0){//reInit to a new direction
                iElm = initElm;
                numElm = elm2ElmConn(iElm, 0);
                if(initDirection == numElm){
                    break;
                }else{
                    initDirection += 1;
                }
                nextEdge = initDirection;
                nextElm = elm2ElmConn(iElm,nextEdge);
                iAcross = -1;
                continue;
            }
            int oldElm = iElm;
            iElm = nextElm;
            numElm = elm2ElmConn(iElm, 0);
            //update the nextEdge:
            for(int i=1; i<=numElm; i++){
                if(elm2ElmConn(iElm, i) == oldElm){
                    nextEdge = i;
                    break;
                }
            }
            nextEdge = (nextEdge + numElm/2)%numElm;
            if(nextEdge == 0)
                nextEdge = numElm;
            nextElm = elm2ElmConn(iElm,nextEdge);
        }
        Vector2 targetPosition = calcElmCenter(iElm, elm2VtxConn, vtxCoords);
        returnDx(iMP) = targetPosition - MPPosition;
    }}); 
    Kokkos::fence();
    
    return returnDx;
}


void runT2LSimple(MPM mpm, int size, double range, const int loopTimes, const int printVTP, const int randomSeed){
    if(size < 0)
        size = mpm.getMPs().getCount();
    if(range < 0)
        range = calcAvgLengthOfEdge(mpm.getMesh())*3;
    printf("Run T2L Tracking Simple Test with size= %d and range= %.4f\n",size, range);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LDelta(size, range, randomSeed);
        mpm.T2LTracking(dx, printVTP<0?-1:i); 
    }
}

void runT2LRankineVortex(MPM mpm, Vector2 center, const int numEdge, double avgLength, const double Gamma, const int loopTimes, const int printVTP){
    if(avgLength <0)
        avgLength = calcAvgLengthOfEdge(mpm.getMesh());
    printf("Run T2L Tracking Rankine Vortex Test with center= (%.4f,%.4f), numEdge= %d, avgLength= %.4f, Gamma= %.4f\n",center[0], center[1], numEdge, avgLength, Gamma);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LDeltaRankineVortex(mpm, center, numEdge, avgLength, Gamma);
        mpm.T2LTracking(dx, printVTP<0?-1:i); 
    }
}

void runT2LRandomWithProportion(MPM mpm, const double p0, const double p1, const double p2, const double p3, const int loopTimes, const int printVTP, const int randomSeed){
    printf("\tRun T2L Tracking Random With Proportion: %.2f %.2f %.2f %.2f\n",p0,p1,p2,p3);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LTest1(mpm, p0, p1, p2, p3, randomSeed);
        mpm.T2LTracking(dx, printVTP<0?-1:i); 
    }
}

void runT2LWithOneMPAcross(MPM mpm, const int MPAcross, const int loopTimes, const int printVTP, const int randomSeed){
    printf("\tRun T2L Tracking with one MP Across %d Element(s)\n",MPAcross);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LTest2(mpm, MPAcross, randomSeed);
        mpm.T2LTracking(dx, printVTP<0?-1:i); 
    }
}

MPM initMPMWithCenterMPs(Mesh& mesh, int factor){
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    auto elm2VtxConn = mesh.getElm2VtxConn(); 

    int numMPs = numElms*factor;
    Vector2View positions("MPpositions", numMPs);
    BoolView isActive("MPstatus",numMPs);
    IntView elm2MPs("elementToMPs",numElms*(maxMPsPerElm+1));
    IntView MPs2Elm("MPToElementIDs",numMPs);

    IntView MPToElement("MPToElement",numMPs);

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
        int iElm = MPToElement(iMP);
        positions(iMP) = calcElmCenter(iElm,elm2VtxConn,vtxCoords);
        MPs2Elm(iMP) = iElm;
        isActive(iMP) = true;
    });    

    auto p = MaterialPoints(numMPs,positions,isActive);
    return MPM(mesh,p,elm2MPs,MPs2Elm);
}

void runCVTRandomWithProportion(MPM mpm, const double p0, const double p1, const double p2, const double p3, const int loopTimes, const int printVTP, const int randomSeed){
    printf("\tRun CVT Tracking Edge Center Based Random With Proportion: %.2f %.2f %.2f %.2f\n",p0,p1,p2,p3);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LTest1(mpm, p0, p1, p2, p3, randomSeed);
        mpm.CVTTrackingEdgeCenterBased(dx, printVTP<0?-1:i); 
    }
}

void runCVTElmCenterBasedRandomWithProportion(MPM mpm, const double p0, const double p1, const double p2, const double p3, const int loopTimes, const int printVTP, const int randomSeed){
    printf("\tRun CVT Tracking Elm Center Based Random With Proportion: %.2f %.2f %.2f %.2f\n",p0,p1,p2,p3);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LTest1(mpm, p0, p1, p2, p3, randomSeed);
        mpm.CVTTrackingElmCenterBased(dx, printVTP<0?-1:i); 
    }
}

void runCVTAllAcrossTest(MPM mpm, const int MPAcross, const double percent, const int loopTimes, const int printVTP, const int randomSeed){
    printf("\tRun CVT Tracking Edge Center Based with all MPs across %d Element(s)\n",MPAcross);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LTest3(mpm, MPAcross, percent, randomSeed);
        mpm.CVTTrackingEdgeCenterBased(dx, printVTP<0?-1:i); 
    }
}

void runCVTElmAllAcrossTest(MPM mpm, const int MPAcross, const double percent, const int loopTimes, const int printVTP, const int randomSeed){
    printf("\tRun CVT Tracking Elm Center Based with all MPs Across %d Element(s)\n",MPAcross);
    for(int i=0; i< loopTimes; i++){
        Vector2View dx = initT2LTest3(mpm, MPAcross, percent, randomSeed);
        mpm.CVTTrackingElmCenterBased(dx, printVTP<0?-1:i); 
    }
}

void printMPs(MPM mpm){
    auto MPs = mpm.getMPs();
    int numMPs = MPs.getCount();
    auto MPsPosition = MPs.getPositions();
    auto MPs2Elm = mpm.getMPs2Elm();
    auto isActive = MPs.isActive();
        
    auto h_MPsPosition = Kokkos::create_mirror_view(MPsPosition);
    auto h_MPs2Elm = Kokkos::create_mirror_view(MPs2Elm);
    auto h_isActive = Kokkos::create_mirror_view(isActive);
    Kokkos::deep_copy(h_MPsPosition,MPsPosition);
    Kokkos::deep_copy(h_MPs2Elm, MPs2Elm);
    Kokkos::deep_copy(h_isActive, isActive);
    printf("The mpm has %d MPs\n",numMPs);
    Kokkos::fence();
    
    for(int i=0; i<numMPs; i++){
        printf("%d: in Elm: %d, at(%.3f,%.3f)\n",i,h_MPs2Elm(i),h_MPsPosition(i)[0],h_MPsPosition(i)[1]);
    }
}
