#include "pmpo_MPMesh.hpp"
#include "pmpo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMPO;

void copyTgtElm(IntView& returnView, MPMesh& mpMesh);

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    const int testMeshOption = 1;
    const int replicateFactor = 1;
    const int testMPOption = 1;
    {//runTracking
        auto testMesh = initTestMesh(testMeshOption, replicateFactor); 
        auto mpMesh = initTestMPMesh(testMesh, testMPOption); 
        auto p_MPs = mpMesh.p_MPs;
        
        Vec2dView dx = polyMPO::Vec2dView("positions",51);
        Kokkos::deep_copy(dx,Vec2d(0.01,0.01));
        
        mpMesh.CVTTrackingEdgeCenterBased(dx);
        IntView TgtElmCVTEdge("TgtElmCVTEdge",p_MPs->getCount());
        copyTgtElm(TgtElmCVTEdge,mpMesh);

        mpMesh.CVTTrackingElmCenterBased(dx);
        IntView TgtElmCVTCenter("TgtElmCVTCenter",p_MPs->getCount());
        copyTgtElm(TgtElmCVTCenter,mpMesh);

        mpMesh.T2LTracking(dx);
        IntView TgtElmT2L("TgtElmT2L",p_MPs->getCount());
        copyTgtElm(TgtElmT2L,mpMesh);

        p_MPs->updateMPSliceAll();
        IntView TgtElmAfter("TgtElmAfter",p_MPs->getCount());
        copyTgtElm(TgtElmAfter,mpMesh);
        //TODO:check all the value
       Kokkos::parallel_for("checkTgtElm",p_MPs->getCount(),KOKKOS_LAMBDA(const int mp){
            //PMT_ALWAYS_ASSERT(TgtElmCVTEdge(mp) == TgtElmCVTCenter(mp)
            //                &&TgtElmCVTEdge(mp) == TgtElmT2L(mp) );
            printf("%d %d %d %d\n",TgtElmCVTEdge(mp),TgtElmCVTCenter(mp),TgtElmT2L(mp),TgtElmAfter(mp));
        });
        
    }
    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

//TODO:change to a template function
void copyTgtElm(IntView& returnView, MPMesh& mpMesh){
    auto mpTgtElm = mpMesh.p_MPs->getData<MPF_Tgt_Elm_ID>();
    auto copyVal = PS_LAMBDA(const int&, const int& mp, const int& mask){
        if(mask){
            returnView(mp) = mpTgtElm(mp);
        }
    };
    mpMesh.p_MPs->parallel_for(copyVal,"copyTgtElmID");
}
