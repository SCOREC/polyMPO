#ifndef PMO_C_H
extern "C" {

typedef void* mpmesh;

void polympo_initialize();
mpmesh polympo_createMpMesh();
void polympo_deleteMpMesh(mpmesh mpMesh);
void polympo_finalize();
void polympo_getMeshCurPosXYZArray(mpmesh mpMeshIn, int size, double* array);
void polympo_setMeshCurPosXYZArray(mpmesh mpMeshIn, int size, double* array);
}
#endif
