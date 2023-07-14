#include "pmt_MPMesh.hpp"
#include "pmt_assembly.hpp"
#include "pmt_wachspressBasis.hpp"

#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMpmTest;

int main(int argc, char* argv[] ) {
    PMT_ALWAYS_ASSERT(argc == 2);
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc,argv);
    int factor = atoi(argv[1]);
    printf("Time assembly and wachspress with factor: %d\n",factor);
    {
        auto mesh = initTestMesh(factor);
        auto mpMesh = initTestMPMesh(mesh,1);
        
        printf("Total MPs:%d\n",mpMesh.MPs->getCount());
        
        for(int i=0; i<5; i++){
            //polyMpmTest::assemblyV0(mpMesh);
            interpolateWachspress(mpMesh);
        }
    } 
    Kokkos::finalize();
    return 0;
}
