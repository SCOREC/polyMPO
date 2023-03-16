#include "MPM.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "assembly.hpp"

#include "testUtils.hpp"

int main(int argc, char* argv[] ) {
    Kokkos::initialize(argc,argv);{
        auto mTest = initTestMesh(atoi(argv[1]));
        auto MPMTest = initTestMPM(mTest);
        
        polyMpmTest::assembly(MPMTest);
        polyMpmTest::assembly(MPMTest);
        polyMpmTest::assembly(MPMTest);
        polyMpmTest::assembly(MPMTest);
        polyMpmTest::assembly(MPMTest);
        
        interpolateWachpress(MPMTest);
        interpolateWachpress(MPMTest);
        interpolateWachpress(MPMTest);
        interpolateWachpress(MPMTest);
        interpolateWachpress(MPMTest);
    } 
    Kokkos::finalize();
    return 0;
}
