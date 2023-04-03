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
        auto mesh = Mesh::readMPASMesh("./grid_island.nc");
        //auto mesh = initTestMesh(1);
        
        auto mpm = initMPMWithRandomMPs(mesh,10);
        
        Vector2View dx = InitT2LDelta(mpm.getMPs().getCount());
        //test T2L in materialPoints
        mpm.T2LTracking(dx);
    }
    
    Kokkos::finalize();
    return 0;
}
