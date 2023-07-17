#ifndef PMO_C_H
#define PMO_C_H

#include <mpi.h>

extern "C" {

typedef void* mpmesh;

//initialize and finalize
void polympo_initialize();
void polympo_finalize();

//create/delete MpMesh object
mpmesh polympo_createMpMesh(int setMeshOption, int setMPOption);
void polympo_deleteMpMesh(mpmesh mpMesh);

//set MPI communicator
void polympo_setMPICommunicator(MPI_Fint fcomm);//TODO:is MPI_Fint best? or something else
//TODO: add a function to get communicator

//MP slices
void polympo_setMPVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_setMP2dVelArray(mpmesh mpMeshIn, int rank1size, int rank2size, double* array);
void polympo_getMPVelArray(mpmesh mpMeshIn, int size, double* array);

//Mesh builds
void polympo_checkMeshSetting(mpmesh mpMeshIn, int maxEdges, int vertexDegree);
void polympo_setMeshNumVtxs(mpmesh mpMeshIn, int numVtxs);
void polympo_setMeshNumElms(mpmesh mpMeshIn, int numElms);
void polympo_setMeshVtxCoords(mpmesh mpMeshIn, int size, double* xArray, double* yArray, double* zArray);
void polympo_setMeshElm2VtxConn(mpmesh mpMeshIn, int size1, int size2, int* array);
void polympo_setMeshVtx2ElmConn(mpmesh mpMeshIn, int size1, int size2, int* array);
void polympo_setMeshElm2ElmConn(mpmesh mpMeshIn, int size1, int size2, int* array);
//Mesh fields
void polympo_setMeshVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_getMeshVelArray(mpmesh mpMeshIn, int size, double* array);
//Mesh helper
void polympo_printMesh(mpmesh mpMeshIn);
}
#endif
