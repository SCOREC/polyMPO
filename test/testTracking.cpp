#include "pmt_MPMesh.hpp"
#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMpmTest;

void copyTgtElm(IntView& returnView, MPMesh& mpMesh);

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    {//runTracking
        auto testMesh = initTestMesh(1); //creates simple test mesh, '1' is a replication factor
        auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
        auto mpMesh = initTestMPMesh(testMesh, mpPerElement); //creates test MPs
        auto MPs = mpMesh.MPs;
        
        Vector2View dx = polyMpmTest::Vector2View("positions",51);
        Kokkos::deep_copy(dx,Vector2(0.01,0.01));
        
        mpMesh.CVTTrackingEdgeCenterBased(dx);
        IntView TgtElmCVTEdge("TgtElmCVTEdge",mpMesh.MPs->getCount());
        copyTgtElm(TgtElmCVTEdge,mpMesh);

        mpMesh.CVTTrackingElmCenterBased(dx);
        IntView TgtElmCVTCenter("TgtElmCVTCenter",mpMesh.MPs->getCount());
        copyTgtElm(TgtElmCVTCenter,mpMesh);

        mpMesh.T2LTracking(dx);
        IntView TgtElmT2L("TgtElmT2L",mpMesh.MPs->getCount());
        copyTgtElm(TgtElmT2L,mpMesh);

        MPs->updateMPSlice();
        IntView TgtElmAfter("TgtElmAfter",mpMesh.MPs->getCount());
        copyTgtElm(TgtElmAfter,mpMesh);
        //TODO:check all the value
       Kokkos::parallel_for("checkTgtElm",MPs->getCount(),KOKKOS_LAMBDA(const int mp){
            //PMT_ALWAYS_ASSERT(TgtElmCVTEdge(mp) == TgtElmCVTCenter(mp)
            //                &&TgtElmCVTEdge(mp) == TgtElmT2L(mp) );
            printf("%d %d %d %d\n",TgtElmCVTEdge(mp),TgtElmCVTCenter(mp),TgtElmT2L(mp),TgtElmAfter(mp));
        });
        
    }
    Kokkos::finalize();
    return 0;
}

//TODO:change to a template function
void copyTgtElm(IntView& returnView, MPMesh& mpMesh){
    auto mpTgtElm = mpMesh.MPs->getData<MPF_Tgt_Elm_ID>();
    auto copyVal = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
        if(mask){
            returnView(mp) = mpTgtElm(mp);
        }
    };
    mpMesh.MPs->parallel_for(copyVal,"copyTgtElmID");
}
