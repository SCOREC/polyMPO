#ifndef PMO_C_H
extern "C" {

typedef void* mpmesh;

void polympo_initialize();
mpmesh polympo_createMpMesh();
void polympo_finalize();
void polympo_modifyArray(int size, double* array);
}
#endif
