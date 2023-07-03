#include "pmt_MPMesh.hpp"
#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMpmTest;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    {//runTracking
        auto testMesh = initTestMesh(1); //creates simple test mesh, '1' is a replication factor
        auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
        auto mpMesh = initTestMPMesh(testMesh, mpPerElement); //creates test MPs
        //auto MPs = mpMesh.MPs;
        
        Vector2View dx = polyMpmTest::Vector2View("positions",51);
        //TODO:set to a value

        mpMesh.CVTTrackingEdgeCenterBased(dx);
        mpMesh.CVTTrackingElmCenterBased(dx);
        mpMesh.T2LTracking(dx);

        //TODO:check all the value
    }
    return 0;
}
