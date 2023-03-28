#include "pmt_MPM.hpp"

#include "testUtils.hpp"

int main() {
    Kokkos::initialize();

    //test T2L
    {
        auto mesh = Mesh::readMPASMesh("./grid_island.nc");
        auto mpm = initTestMPM(mesh);
        
        Vector2View dx = InitT2LDelta(mpm.getMPs().getCount());
        //test T2L in materialPoints
        mpm.T2LTracking(dx);
    }
    
    Kokkos::finalize();
    return 0;
}