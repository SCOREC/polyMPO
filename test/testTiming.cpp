#include "MPM.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "assembly.hpp"

#include "testUtils.hpp"

int main(int argc, char* argv[] ) {
    int factor = atoi(argv[1]);
    printf("Time assembly and wachspress with factor: %d\n",factor);
    Kokkos::initialize(argc,argv);{
        auto mTest = initTestMesh(factor);
        auto MPMTest = initTestMPM(mTest);
        
        printf("Total MPs:%d\n",MPMTest.getMPs().getCount());    
    
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
