#ifndef POLYMPMTEST_ASSEMBLY_H
#define POLYMPMTEST_ASSEMBLY_H

namespace polyMpmTest{

DoubleView assembly(MPM& mpm){
    auto mesh = mpm.getMesh();
    auto MPs = mpm.getMPs();
    int numVtxs = mesh.getNumVertices();
    int numElms = mesh.getNumElements();
    int numMPs = MPs.getCount(); 
     
    auto vtxCoords = mesh.getVtxCoords(); 
    auto elm2VtxConn = mesh.getElm2VtxConn();
    auto xp = MPs.getPositions();
    
    auto MPs2Elm = mpm.getMPs2Elm();
    DoubleView vField("vField2",numVtxs);
    Kokkos::parallel_for("vertex_assem", numMPs, KOKKOS_LAMBDA(const int iMP){
        int ielm = MPs2Elm(iMP); 
        int nVtxE = elm2VtxConn(ielm,0);
        for(int i=0; i<nVtxE; i++){
            int vID = elm2VtxConn(ielm,i+1)-1;
            auto vertexLoc = vtxCoords(vID);
            double distance = (xp(iMP)-vertexLoc).magnitude();
            Kokkos::atomic_add(&vField(vID),distance);
        }
    });
   
    //Kokkos::parallel_for("vFieldcheck",numVtxs, KOKKOS_LAMBDA(const int ielm){
    //    printf("%d: %f %f\n",ielm,vField(ielm),vField2(ielm));
    //});
    return vField;
}

} //end namespace polyMpmTest
#endif
