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

Vector2View initT2LTest1(MPM mpm, double percent1, double percent2, double percent3, double percent4, const int randomSeed){
    PMT_ALWAYS_ASSERT(percent1+percent2+percent3+percent4 <= 1+MPMTEST_EPSILON);
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();
    Vector2View vtxCoords = mesh.getVtxCoords();   
    int numMPs = MPs.getCount();
    auto MPsPosition = MPs.getPositions();
    auto MPs2Elm = mpm.getMPs2Elm();
    auto elm2ElmConn = mesh.getElm2ElmConn();
    auto isActive = MPs.isActive();
    auto elm2VtxConn = mesh.getElm2VtxConn();
    percent2 += percent1;
    percent3 += percent2;
    percent4 += percent3;
    Vector2View returnDx("T2LDeltaX",numMPs);
    Kokkos::Random_XorShift64_Pool<> random_pool(randomSeed);
    Kokkos::parallel_for("setNumMPPerElement", numMPs, KOKKOS_LAMBDA(const int iMP){
    if(isActive(iMP)){
        auto generator = random_pool.get_state();
        Vector2 MPPosition = MPsPosition(iMP);
        double iRange = generator.drand(1);
        int numAcross = 0;
        if(iRange < percent1){
        }else if(iRange < percent2){
            numAcross = 1;
        }else if(iRange < percent3){
            numAcross = 2;
        }else if(iRange < percent4){
            numAcross = 3;
        }else{
            numAcross = 4;
        }
        int initElm = MPs2Elm(iMP);
        int iElm = initElm;
        Vector2 currentCenter = calcElmCenter(iElm,elm2VtxConn,vtxCoords);
        int numEdge = elm2ElmConn(iElm, 0);
        int numVtx = elm2VtxConn(iElm,0);
    // calc the center line:
        int nextEdge = 0;
        int nextElm = 0;
        Vector2 MP = MPsPosition(iMP);
        Vector2 dx = (currentCenter - MP)*100000;//arbitrary size only get the direction
        int v[maxVtxsPerElm];
        for(int i=0; i< numVtx; i++)
            v[i] = elm2VtxConn(iElm,i+1)-1;
        Vector2 e[maxVtxsPerElm];
        double pdx[maxVtxsPerElm];                    
        for(int i=0; i< numVtx; i++){
            Vector2 v_i = vtxCoords(v[i]);
            Vector2 v_ip1 = vtxCoords(v[(i+1)%numVtx]);
            e[i] = v_ip1 - v_i;
            pdx[i] = (v_i - MP).cross(dx);
        }
        for(int i=0; i<numVtx; i++){
            int ip1 = (i+1)%numVtx;
        //printf("%d,%d (%f,%f)\n",nextElm,nextEdge,dx[0],dx[1]);
            if(pdx[i]*pdx[ip1] <0){
            //iElm = elm2ElmConn(iElm,i+1);
                nextElm = elm2ElmConn(iElm,i+1);
                for(int j=1; j<=numEdge; j++){
                    if(elm2ElmConn(iElm, j) == initElm){
                        nextEdge = j;
                        break;
                    }
                }
            //goToNeighbour = true;
            break;
            }
        }
    //
        //int initDirection = generator.urand(0,numEdge)+1;//(0:numEdge)+1
        int initDirection =  nextEdge;
        //int revDirection = (initDirection + numEdge/2)%numEdge;
        //bool goOut = false;
        //bool goRev = false;
        //double printfx[2] = {563196,563197};
        for(int iAcross = 0; iAcross< numAcross; iAcross++){
            //if(MPPosition[0] <printfx[1] && MPPosition[0] > printfx[0]){
            //    printf("i(%d): (%f,%f), initElm: %d, iElm: %d, nextElm: %d\n",numAcross,vtxCoords(elm2VtxConn(iElm,1)-1)[0],vtxCoords(elm2VtxConn(iElm,1)-1)[1],initElm,iElm,nextElm);
            //}
            //if(MPPosition[0] <-629065 && MPPosition[0] >-629066){
            //    printf("%d:nextEdge: %d, nextElm: %d\n",iMP,nextEdge,nextElm);
            //}
            /*if(nextElm < 0){
                if(!goRev){
                    goRev = true;
                    nextEdge = revDirection;
                    nextElm = elm2ElmConn(iElm,nextEdge);
                    if(nextElm < 0){
                        goOut = true;
                        break;
                    }
                    iAcross = -1;
                    continue;
                }else{
                    goOut = true;
                    //printf("go out\n");
                    break;
                }
            }*/
            if(nextElm < 0){//reInit to a new direction
                iElm = initElm;
                numEdge = elm2ElmConn(iElm, 0);
                if(initDirection == numEdge){
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
            numEdge = elm2ElmConn(iElm, 0);
            //update the nextEdge:
            for(int i=1; i<=numEdge; i++){
                if(elm2ElmConn(iElm, i) == oldElm){
                    nextEdge = i;
                    break;
                }
            }
            //if(MPPosition[0] <printfx[1] && MPPosition[0] > printfx[0]){
            //    printf("%d -> ",nextEdge);
            //}
            nextEdge = (nextEdge + numEdge/2)%numEdge;
            //if(MPPosition[0] <printfx[1] && MPPosition[0] > printfx[0]){
            //    printf("%d\n",nextEdge);
            //}
            if(nextEdge == 0)
                nextEdge = numEdge;
            nextElm = elm2ElmConn(iElm,nextEdge);
        }
            //if(MPPosition[0] <printfx[1] && MPPosition[0] > printfx[0]){
            //    printf("i(%d): (%f,%f), initElm: %d, iElm: %d, nextElm: %d\n",numAcross,vtxCoords(elm2VtxConn(iElm,1)-1)[0],vtxCoords(elm2VtxConn(iElm,1)-1)[1],initElm,iElm,nextElm);
            //}
        //normal case nextElm >= 0;
        Vector2 XYc = calcElmCenter(iElm,elm2VtxConn, vtxCoords);
        //int triID = generator.urand(0,numVtx);
        int triID = nextEdge-1;
        //printf("nextElm= %d, nextEdge= %d, nextElmFromNextEdge= %d\n",nextElm, nextEdge, elm2ElmConn(iElm,nextEdge));
        double rws[2] = {generator.drand(0.0,1.0), generator.drand(0.0,1.0)};
        random_pool.free_state(generator);
        //double weights[3];
        //if (rws[0]> rws[1]){
        //    weights[0] = rws[1];
        //    weights[1] = rws[0]-rws[1];
        //    weights[2] = 1-rws[0];
        //}else{
        //    weights[0] = rws[0];
        //    weights[1] = rws[1]-rws[0];
        //    weights[2] = 1-rws[1];
        //}
        auto v1 = vtxCoords(elm2VtxConn(iElm,triID+1)-1);
        auto v2 = vtxCoords(elm2VtxConn(iElm,(triID+1)%numVtx+1)-1);
        //Vector2 targetPosition = XYc*weights[0]+v1*weights[1]+v2*weights[2];
        Vector2 targetPosition = XYc;
        //go out nextElm <0;
        /*if(goOut){
            //find the two vertex then -y, x
            if(nextEdge == elm2VtxConn(iElm,0)){
                v1 = vtxCoords(elm2VtxConn(iElm, nextEdge)-1);
                v2 = vtxCoords(elm2VtxConn(iElm, 1)-1);
            }else{
                v1 = vtxCoords(elm2VtxConn(iElm, nextEdge)-1);
                v2 = vtxCoords(elm2VtxConn(iElm, nextEdge+1)-1);
            }
            Vector2 diff = v1-v2;
            targetPosition = targetPosition + Vector2(diff[1], -diff[0]);
            //targetPosition = targetPosition*1.1;
        }*/     
        //if(MPPosition[0] <printfx[1] && MPPosition[0] > printfx[0]){
        //    printf("%d: initElm: %d, iElm: %d, v1= (%f,%f), v2= (%f,%f)\n",goOut,initElm,iElm,v1[0],v1[1],v2[0],v2[1]);
        //}
        returnDx(iMP) = targetPosition - MPPosition;
        //if(iMP < 10)
        //    printf("%d:(%f) = (%f,%f)-(%f,%f)\n", iMP, returnDx(iMP).magnitude(), targetPosition[1], targetPosition[1], MPPosition[0], MPPosition[1]);
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
