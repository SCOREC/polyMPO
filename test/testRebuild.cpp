#include "pmt_MPMesh.hpp"
#include "pmo_createTestMPMesh.hpp"
#include "testUtils.hpp"
#include <mpi.h>

using namespace polyMpmTest;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
    auto testMesh = initTestMesh(1); //creates simple test mesh, '1' is a replication factor
    auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
    auto mpMesh = initTestMPMesh(testMesh, mpPerElement); //creates test MPs
    auto MPs = mpMesh.MPs;

    //move the mps to their current element id % 2
    auto mpTgtElm = MPs->getData<MPF_Tgt_Elm_ID>();
    auto mpStatus = MPs->getData<MPF_Status>();
    auto setTargetElmMod2 = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) { 
        mpStatus(mp) = elm;
        mpTgtElm(mp) = elm % 2;
      }
    };
    MPs->parallel_for(setTargetElmMod2, "setTargetElmMod2");
    MPs->rebuild();
    auto checkTargetElmMod2 = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) assert(mpTgtElm(mp) == mpStatus(mp) % 2);
    };
    MPs->parallel_for(checkTargetElmMod2, "setTargetElmMod2");

    //move the all mps to elmeent zero, the zeroth mp is marked as inactive
    const auto initNumMPs = MPs->getCount();
    auto setTargetElmZero = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) { 
        mpTgtElm(mp) = (mp==0) ? -1 : 0;
      }
    };
    MPs->parallel_for(setTargetElmZero, "setTargetElmZero");
    MPs->rebuild();
    PMT_ALWAYS_ASSERT(MPs->getCount() == initNumMPs - 1);
    auto checkTargetElmZero = PS_LAMBDA(const int& elm, const int& mp, const int& mask){
      if(mask) assert(mpTgtElm(mp) == 0);
    };
    MPs->parallel_for(checkTargetElmZero, "setTargetElmZero");

  }
  Kokkos::finalize();

  return 0;
}
