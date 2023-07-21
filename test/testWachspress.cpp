#include "pmt_MPMesh.hpp"
#include "pmt_wachspressBasis.hpp"
#include "pmt_assembly.hpp"

#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"

using namespace polyMpmTest;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    const int testMeshOption = 1;
    const int scaleFactor = 100;
    const int testMPOption = 1;
    //test init Test Mesh and run assembly and Wachspress
    {
        auto m = polyMpmTest::Mesh();
        auto mp = polyMpmTest::MaterialPoints();
        auto v = polyMpmTest::Vec2d();
        
        auto mesh = initTestMesh(testMeshOption, scaleFactor);
        auto mpMesh = initTestMPMesh(mesh,testMPOption);
        
        //test assembly in assembly.hpp
        polyMpmTest::assembly<MPF_Vel,MeshF_Vel>(mpMesh,false,false);
        interpolateWachspress(mpMesh);
    }
    
    Kokkos::finalize();
    return 0;
}
