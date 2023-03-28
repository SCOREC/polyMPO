#include "pmt_MPM.hpp"
#include "pmt_wachspressBasis.hpp"
#include "pmt_assembly.hpp"

#include "testUtils.hpp"

int main() {
    Kokkos::initialize();
    
    //test init Test Mesh and run assembly and Wachspress
    {
        auto m = polyMpmTest::Mesh();
        auto mp = polyMpmTest::MaterialPoints();
        auto v = polyMpmTest::Vector2();
        
        auto mesh = initTestMesh(100);
        auto mpm = initTestMPM(mesh);
        
        //test assembly in assembly.hpp
        polyMpmTest::assembly(mpm);
        interpolateWachspress(mpm);
    }
    
    Kokkos::finalize();
    return 0;
}
