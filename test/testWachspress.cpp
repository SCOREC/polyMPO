#include "pmpo_MPMesh.hpp"
#include "pmpo_wachspressBasis.hpp"
#include "pmpo_MPMesh_assembly.hpp"
#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMPO;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    const int replicateFactor = 100;
    const int testMPOption = 1;
    //test init Test Mesh and run assembly and Wachspress
    for (int testMeshOption = 1; testMeshOption <= 2; testMeshOption++){

        auto mesh = initTestMesh(testMeshOption,replicateFactor);
        auto mpMesh = initTestMPMesh(mesh,testMPOption);
        
        //test assembly in assembly.hpp
        mpMesh.assembly<MeshF_Vel>(0, polyMPO::MeshFType_VtxBased, false,false);
        interpolateWachspress2DTest(mpMesh);
        interpolateWachspress3DTest(mpMesh,testMeshOption);
    }
    
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}
