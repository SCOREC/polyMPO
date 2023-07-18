#include "pmt_MPMesh.hpp"
#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMpmTest;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
 
  const int testMeshOption = 1;
  const int scaleFactor = 1;
  const int testMPOption = 1;
  {
    polyMpmTest::Mesh* testMesh = initTestMesh(testMeshOption, scaleFactor); 
    auto mpMesh = initTestMPMesh(testMesh, testMPOption); //creates test MPs
    auto p_MPs = mpMesh.p_MPs;

    //move the mps to their current element id % 2
    auto mpTgtElm = p_MPs->getData<MPF_Tgt_Elm_ID>();
    auto mpStatus = p_MPs->getData<MPF_Status>();
    auto setTargetElmMod2 = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) { 
        mpStatus(mp) = elm;
        mpTgtElm(mp) = elm % 2;
      }
    };
    p_MPs->parallel_for(setTargetElmMod2, "setTargetElmMod2");
    p_MPs->rebuild();
    auto checkTargetElmMod2 = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) assert(mpTgtElm(mp) == mpStatus(mp) % 2);
    };
    p_MPs->parallel_for(checkTargetElmMod2, "setTargetElmMod2");

    //move the all mps to elmeent zero, the zeroth mp is marked as inactive
    const auto initNumMPs = p_MPs->getCount();
    auto setTargetElmZero = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) { 
        mpTgtElm(mp) = (mp==0) ? -1 : 0;
      }
    };
    p_MPs->parallel_for(setTargetElmZero, "setTargetElmZero");
    p_MPs->rebuild();
    PMT_ALWAYS_ASSERT(p_MPs->getCount() == initNumMPs - 1);
    auto checkTargetElmZero = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) assert(mpTgtElm(mp) == 0);
    };
    p_MPs->parallel_for(checkTargetElmZero, "setTargetElmZero");

  }
  Kokkos::finalize();

  return 0;
}
