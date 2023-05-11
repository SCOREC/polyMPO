#include "pmt_MPMesh.hpp"
#include "pmt_utils.hpp"
#include <mpi.h>

void checkPositions(polyMpmTest::MaterialPoints& MPs, std::string name) {
  auto mpPositions = MPs.getData<0>();
  auto checkPositions = PS_LAMBDA(const int&, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      PMT_ALWAYS_ASSERT(mpPositions(mp,0) == 1.0);
      PMT_ALWAYS_ASSERT(mpPositions(mp,1) == 42.0);
    }
  };
  MPs.parallel_for(checkPositions, "checkPositions_"+name);
}


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
    const int numElms = 10;
    const int numMPs = 20;
    auto mpPositions = polyMpmTest::Vector2View("positions",numMPs);
    Kokkos::parallel_for("intializeMPsPosition", numMPs, KOKKOS_LAMBDA(const int i){
        const auto x = 1.0;
        const auto y = 42.0;
        mpPositions(i) = polyMpmTest::Vector2(x,y);
    });
    auto mpPerElement = std::vector<int>({2,3,1,2,3,3,2,1,2,1});
    const auto mpCount = std::accumulate(mpPerElement.begin(),mpPerElement.end(),0);
    PMT_ALWAYS_ASSERT(mpCount == numMPs);
    auto MPs = polyMpmTest::MaterialPoints(numElms, numMPs, mpPositions);
    checkPositions(MPs, "afterConstruction");

    auto mpToElement = std::vector<int>(numMPs);
    int mpIdx = 0;
    for(int i=0; i<numElms; i++) {
      for(int j=0; j<mpPerElement[i]; j++) {
        mpToElement[mpIdx++] = i;
      }
    }
    Kokkos::View<int*, Kokkos::HostSpace> mpToElement_h(mpToElement.data(),mpToElement.size());
    polyMpmTest::IntView mpToElement_d("mpToElement",numMPs);
    Kokkos::deep_copy(mpToElement_d,mpToElement_h);
    MPs.rebuild(mpToElement_d);
    checkPositions(MPs, "afterRebuild");
  }
  Kokkos::finalize();
  return 0;
}
