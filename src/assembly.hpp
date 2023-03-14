#ifndef POLYMPMTEST_ASSEMBLY_H
#define POLYMPMTEST_ASSEMBLY_H

namespace polyMpmTest{

void assembly(MPM mpm){
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    int numMPs = MPs.getCount(); 
     
    auto vtxCoords = mesh.getVtxCoords(); 
    auto elm2VtxConn = mesh.getElm2VtxConn();
    auto elm2MPs = mpm.getElm2MPs();
    auto xp = MPs.getPositions();

    DoubleView vField("vField",numVtxs);
    
    Kokkos::parallel_for("vertex_assem",numElms, KOKKOS_LAMBDA(const int ielm){
        int nVtxE = elm2VtxConn(ielm,0);
        int nMPE = elm2MPs(ielm*(maxMPsPerElm+1)); 
        for(int i=0; i<nVtxE; i++){
            int vID = elm2VtxConn(ielm,i+1)-1;
            auto vertexLoc = vtxCoords(vID);
            for(int j=0; j<nMPE; j++){
                int MPID = elm2MPs(ielm*(maxMPsPerElm+1)+j+1);
                double distance = (xp(MPID)-vertexLoc).magnitude();
                Kokkos::atomic_add(&vField(vID),distance);
            }
        }
    });
    auto MPs2Elm = mpm.getMPs2Elm();
    DoubleView vField2("vField2",numVtxs);
    Kokkos::parallel_for("vertex_assem2", numMPs, KOKKOS_LAMBDA(const int iMP){
        int ielm = MPs2Elm(iMP); 
        int nVtxE = elm2VtxConn(ielm,0);
        for(int i=0; i<nVtxE; i++){
            int vID = elm2VtxConn(ielm,i+1)-1;
            auto vertexLoc = vtxCoords(vID);
            double distance = (xp(iMP)-vertexLoc).magnitude();
            Kokkos::atomic_add(&vField2(vID),distance);
        }
    });
   
    //Kokkos::parallel_for("vFieldcheck",numVtxs, KOKKOS_LAMBDA(const int ielm){
    //    printf("%d: %f %f\n",ielm,vField(ielm),vField2(ielm));
    //});

}

} //end namespace polyMpmTest
#endif
