#include "pmo_c.h"
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
void polympo_modifyArray(int size, double* array) {
  printf("polympo_modifyArray c++ size %d\n", size);
  array[0] = 42.0;
}

#ifdef __cplusplus
}
#endif
