#include "pmo_c.h"
#include "pmo_createTestMPMesh.hpp"
#include <stdio.h>


void polympo_initialize() {
  printf("polympo_initialize c++\n");
  int isMpiInit;
  MPI_Initialized(&isMpiInit);
  PMT_ALWAYS_ASSERT(isMpiInit);
  Kokkos::initialize();
}

void polympo_finalize() {
  printf("polympo_finalize c++\n");
  Kokkos::finalize();
}

mpmesh polympo_createMpMesh() {
  auto mesh = polyMpmTest::initTestMesh(1);
  auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
  auto mps = polyMpmTest::initTestMPs(mesh, mpPerElement);
  return (mpmesh) new polyMpmTest::MPMesh(mesh,mps);
}

void polympo_deleteMpMesh(mpmesh mpMeshIn) {
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  delete mpMesh;
}

/**
 * Attention: this typedef is LayoutLeft, meaning that the first
 * index is the contiguous one. This matches the Fortran and GPU conventions for
 * allocations.
 */
typedef Kokkos::View<
          double*[3],
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>
        > kkDblViewHostU;

void polympo_setMeshCurPosXYZArray(mpmesh mpMeshIn, int size, double* array) {
  printf("polympo_setMeshArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = mpMesh->getMesh();
  auto vtxField = mesh.getMeshField<polyMpmTest::MeshF_Cur_Pos_XYZ>();

  //check the size
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*3)==vtxField.size());

  //copy the host array to the device
  Kokkos::deep_copy(vtxField,arrayHost);

  //modify the array - this is just for testing
  Kokkos::parallel_for("setOneEntry", 1, KOKKOS_LAMBDA(int) {
      vtxField(0,0) = 42.0;
  });
}

void polympo_getMeshCurPosXYZArray(mpmesh mpMeshIn, int size, double* array) {
  printf("polympo_getMeshArray c++ size %d\n", size);
  polyMpmTest::MPMesh* mpMesh = (polyMpmTest::MPMesh*)mpMeshIn;
  kkDblViewHostU arrayHost(array,size);

  auto mesh = mpMesh->getMesh();
  auto vtxField = mesh.getMeshField<polyMpmTest::MeshF_Cur_Pos_XYZ>();

  //check the size
  PMT_ALWAYS_ASSERT(static_cast<size_t>(size*3)==vtxField.size());

  Kokkos::deep_copy(arrayHost, vtxField);
}
