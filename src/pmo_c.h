#ifndef PMO_C_H
#define PMO_C_H

#include <mpi.h>

extern "C" {

typedef void* mpmesh;

void polympo_initialize();
mpmesh polympo_createMpMesh();
void polympo_deleteMpMesh(mpmesh mpMesh);
void polympo_finalize();

void polympo_setMPVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_setMP2dVelArray(mpmesh mpMeshIn, int rank1size, int rank2size, double* array);
void polympo_getMPVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_setCommunicator(MPI_Fint fcomm);

void polympo_setMeshVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_getMeshVelArray(mpmesh mpMeshIn, int size, double* array);

void polympo_setMeshVtxCoords(mpmesh mpMeshIn, int size, double* xArray, double* yArray, double* zArray);
void polympo_setMeshElm2VtxConn(mpmesh mpMeshIn, int size, int* array);
void polympo_setMeshVtx2ElmConn(mpmesh mpMeshIn, int size, int* array);
void polympo_setMeshElm2ElmConn(mpmesh mpMeshIn, int size1, int size2, int* array);

void polympo_checkMeshSetting(mpmesh mpMeshIn,int maxEdges,int vertexDegree);
void polympo_setNumVtxs(mpmesh mpMeshIn, int numVtxs);
void polympo_setNumElms(mpmesh mpMeshIn, int numElms);
void polympo_printMesh(mpmesh mpMeshIn);
}
#endif
