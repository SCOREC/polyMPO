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
        Kokkos::deep_copy(dx,Vector2(100,100));
        
        mpMesh.CVTTrackingEdgeCenterBased(dx);
        mpMesh.CVTTrackingElmCenterBased(dx);
        mpMesh.T2LTracking(dx);

        //TODO:check all the value
        auto checkVal = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
            //if(mask){
                printf("mp:%d in Elm:%d is %d\n",mp,elm,mask);
            //}
        };
        mpMesh.MPs->parallel_for(checkVal,"test Tracking check values");
    }
    return 0;
}
