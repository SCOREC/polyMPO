#include "pmt_MPMesh.hpp"
#include "pmt_wachspressBasis.hpp"
#include "pmt_assembly.hpp"

#include "testUtils.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    
    //test init Test Mesh and run assembly and Wachspress
    {
        auto m = polyMpmTest::Mesh();
        auto mp = polyMpmTest::MaterialPoints();
        auto v = polyMpmTest::Vector2();
        
        auto mesh = initTestMesh(100);
        auto mpMesh = initTestMPMesh(mesh);
        
        //test assembly in assembly.hpp
        polyMpmTest::assembly(mpMesh);
        interpolateWachspress(mpMesh);
    }
    
    Kokkos::finalize();
    return 0;
}
