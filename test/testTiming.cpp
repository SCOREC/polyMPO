#include "pmpo_MPMesh.hpp"
#include "pmpo_assembly.hpp"
#include "pmpo_wachspressBasis.hpp"

#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

int main(int argc, char* argv[] ) {
    PMT_ALWAYS_ASSERT(argc == 2);
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc,argv);
    const int testMeshOption = 1;
    const int scaleFactor = atoi(argv[1]);
    const int testMPOption = 1;
    printf("Time assembly and wachspress with factor: %d\n",scaleFactor);
    {
        auto mesh = initTestMesh(testMeshOption,scaleFactor);
        auto mpMesh = initTestMPMesh(mesh,testMPOption);
        
        printf("Total MPs:%d\n",mpMesh.p_MPs->getCount());
        
        for(int i=0; i<5; i++){
            //polyMPO::assemblyV0(mpMesh);
            interpolateWachspress(mpMesh);
        }
    } 
    Kokkos::finalize();
    return 0;
}
