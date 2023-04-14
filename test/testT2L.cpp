#include "pmt_MPM.hpp"

#include "testUtils.hpp"

int main() {
    Kokkos::initialize();
    
    //for(int i=0; i<8704; i++){
    //    printf("%d %d\n",i*2, i*2+1);
    //}
    //test T2L
    {
        //auto mesh = Mesh::readMPASMesh("./mesh.QU.1920km.151026.nc");
        auto mesh = Mesh::readMPASMesh("./GIS.nc");
        //auto mesh = initTestMesh(1);
        
        auto mpm = initMPMWithRandomMPs(mesh,10);
        
        //Vector2View dx = initT2LDelta(mpm.getMPs().getCount(),4000,12345);
        //test T2L in materialPoints
        //mpm.T2LTracking(dx);

        initT2LTest2(mpm.getMPs().getCount(),11500,0.7,0.1,0.1,0.1);
        //calcAvgLengthOfEdge(mesh,1); 
        //for(int i=0; i<100; i++){ 
            //TODO: improve the writing to different file
        //    Vector2View dx = initT2LDeltaRankineVortex(mpm, Vector2(150000, -2000000), 15, 11500, 1);
        //    mpm.T2LTracking(dx,-1);
        //}
    }    
    Kokkos::finalize();
    return 0;
}
