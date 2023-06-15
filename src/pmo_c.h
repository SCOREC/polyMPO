#ifndef PMO_C_H
extern "C" {

typedef void* mpmesh;

void polympo_initialize();
mpmesh polympo_createMpMesh();
void polympo_deleteMpMesh(mpmesh mpMesh);
void polympo_finalize();
void polympo_setMPVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_getMPVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_setMeshVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_getMeshVelArray(mpmesh mpMeshIn, int size, double* array);
}
#endif
