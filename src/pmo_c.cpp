#include "pmo_c.h"
#include "pmo_createTestMPMesh.hpp"
#include <stdio.h>

void polympo_initialize() {
  printf("polympo_initialize c++\n");
}
void polympo_finalize() {
  printf("polympo_finalize c++\n");
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

void polympo_modifyArray(int size, double* array) {
  printf("polympo_modifyArray c++ size %d\n", size);
  array[0] = 42.0;
}
