#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"
#include "pmpo_assembly.hpp"
#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

void test_planar(const int replicateFactor, const int testMPOption) {
    //test init Test Mesh and run assembly and Wachspress
    for (int testMeshOption = 1; testMeshOption <= 2; testMeshOption++){

        auto mesh = initTestMesh(testMeshOption,replicateFactor);
        auto mpMesh = initTestMPMesh(mesh,testMPOption);
        
        //test assembly in assembly.hpp
        polyMPO::assembly<MPF_Vel,MeshF_Vel>(mpMesh,false,false);
        interpolateWachspress2DTest(mpMesh);
        interpolateWachspress3DTest(mpMesh);
        //interpolateWachspressSphericalTest(mpMesh);
        interpolateWachspress3DTest(mpMesh,testMeshOption);
    }

}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    const int replicateFactor = 100;
    const int testMPOption = 1;
        
    test_planar(replicateFactor, testMPOption);
    
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
