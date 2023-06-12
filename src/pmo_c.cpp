#include "pmo_c.h"
#include "testUtils.hpp"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void polympo_initialize() {
  printf("polympo_initialize c++\n");
}
void polympo_finalize() {
  printf("polympo_finalize c++\n");
}

mpmesh polympo_createMpMesh() {
  auto testMesh = initTestMesh(1); //creates simple test mesh, '1' is a replication factor
  auto mpPerElement = std::vector<int>({5,4,5,6,6,5,4,6,5,5});
  auto mpMesh = initTestMPMesh(testMesh, mpPerElement); //creates test MPs
  return (mpmesh)mpMesh;//TODO - need to return a pointer - not going to define custom fortran type....
}
void polympo_modifyArray(int size, double* array) {
  printf("polympo_modifyArray c++ size %d\n", size);
  array[0] = 42.0;
}

#ifdef __cplusplus
}
#endif
