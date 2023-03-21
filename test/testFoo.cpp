#include "MPM.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "assembly.hpp"

#include "testUtils.hpp"

int main() {
    Kokkos::initialize();

    //test init
    auto m = polyMpmTest::Mesh();
    auto mp = polyMpmTest::MaterialPoints();
    auto v = polyMpmTest::Vector2();
    
    //test init Test Mesh and run assembly and Wachspress
    {
        auto mesh = initTestMesh(100);
        auto mpm = initTestMPM(mesh);
        
        //test assembly in assembly.hpp
        polyMpmTest::assembly(mpm);
        interpolateWachpress(mpm);
    }
    
    Kokkos::finalize();
    return 0;
}
