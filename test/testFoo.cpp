#include "MPM.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "assembly.hpp"

#include "testUtils.hpp"

int main() {

    //test init
    auto m = polyMpmTest::Mesh();
    auto mp = polyMpmTest::MaterialPoints();
    auto v = polyMpmTest::Vector2();
    
    //test init Test Mesh
    Kokkos::initialize();{
        auto mTest = initTestMesh(100);
        auto MPMTest = initTestMPM(mTest);
        
        
            
        //test assembly in assembly.hpp
        polyMpmTest::assembly(MPMTest);
        interpolateWachpress(MPMTest);
    }
    
    Kokkos::finalize();
    return 0;
}
