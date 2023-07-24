#include "pmt_MPMesh.hpp"
#include "pmt_utils.hpp"
#include <mpi.h>

void checkPositions(polyMpmTest::MaterialPoints& MPs, std::string name) {
  auto mpPositions = MPs.getData<polyMpmTest::MPF_Cur_Pos_XYZ>();
  auto checkPositions = PS_LAMBDA(const int& elm, const int& mp, const int& mask) {
    if(mask) { //if material point is 'active'/'enabled'
      assert(mpPositions(mp,0) == elm*1.0);
      assert(mpPositions(mp,1) == elm*-1.0);
    }
  };
  MPs.parallel_for(checkPositions, "checkPositions_"+name);
}

template <typename T>
auto copyToDevice(std::vector<T> inVec, std::string name) {
  Kokkos::View<int*, Kokkos::HostSpace> src(inVec.data(),inVec.size());
  Kokkos::View<T*, Kokkos::DefaultExecutionSpace> dest(name,inVec.size());
  Kokkos::deep_copy(dest,src);
  return dest;
}


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
    const int numElms = 10;
    const int numMPs = 20;

    auto mpPerElement = std::vector<int>({2,3,1,2,3,3,2,1,2,1});
    const auto mpCount = std::accumulate(mpPerElement.begin(),mpPerElement.end(),0);
    PMT_ALWAYS_ASSERT(mpCount == numMPs);

    auto mpToElement = std::vector<int>(numMPs);
    int mpIdx = 0;
    for(int i=0; i<numElms; i++) {
      for(int j=0; j<mpPerElement[i]; j++) {
        mpToElement[mpIdx++] = i;
      }
    }

    auto mpToElement_d = copyToDevice<int>(mpToElement, "mpToElement");
    auto mpPerElement_d = copyToDevice<int>(mpPerElement, "mpPerElement");

    auto mpPositions = polyMpmTest::DoubleVec3dView("positions",numMPs);
    Kokkos::parallel_for("intializeMPsPosition", numMPs, KOKKOS_LAMBDA(const int i){
        const auto x = mpToElement_d(i)*1.0;
        const auto y = mpToElement_d(i)*-1.0;
        mpPositions(i,0) = x;
        mpPositions(i,1) = y;
        mpPositions(i,2) = 0.0;
    });

    auto MPs = polyMpmTest::MaterialPoints(numElms, numMPs, mpPositions, mpPerElement_d, mpToElement_d);
    checkPositions(MPs, "afterConstruction");
  }
  Kokkos::finalize();
  return 0;
}
