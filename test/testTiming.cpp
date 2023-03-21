#include "MPM.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "assembly.hpp"

#include "testUtils.hpp"

int main(int argc, char* argv[] ) {
    PMT_ALWAYS_ASSERT(argc == 1);
    Kokkos::initialize(argc,argv);
    int factor = atoi(argv[1]);
    printf("Time assembly and wachspress with factor: %d\n",factor);
    {
        auto mesh = initTestMesh(factor);
        auto mpm = initTestMPM(mesh);
        
        printf("Total MPs:%d\n",mpm.getMPs().getCount());    
        
        for(int i=0; i<5; i++){
            polyMpmTest::assembly(mpm);
            interpolateWachpress(mpm);
        }
    } 
    Kokkos::finalize();
    return 0;
}
