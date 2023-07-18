#ifndef PMO_C_H
#define PMO_C_H

#include <mpi.h>

extern "C" {

typedef void* MPMesh_ptr;

//initialize and finalize
void polympo_initialize();
void polympo_finalize();

//create/delete MpMesh object
MPMesh_ptr polympo_createMPMesh(int setMeshOption, int setMPOption);
void polympo_deleteMPMesh(MPMesh_ptr p_mpmesh);

//set MPI communicator
void polympo_setMPICommunicator(MPI_Fint fcomm);//TODO:is MPI_Fint best? or something else
//TODO: add a function to get communicator

//MP info
//TODO 

//MP slices
void polympo_setMPVelArray(MPMesh_ptr p_mpmesh, int size, double* array);
void polympo_getMPVelArray(MPMesh_ptr p_mpmesh, int size, double* array);

//Mesh info
void polympo_checkMeshSetting(MPMesh_ptr p_mpmesh, int maxEdges, int vertexDegree);
void polympo_setMeshNumVtxs(MPMesh_ptr p_mpmesh, int numVtxs);
void polympo_setMeshNumElms(MPMesh_ptr p_mpmesh, int numElms);
void polympo_setMeshVtxCoords(MPMesh_ptr p_mpmesh, int size, double* xArray, double* yArray, double* zArray);
void polympo_setMeshNumEdgesPerElm(MPMesh_ptr p_mpmesh, int size, int* array);
void polympo_setMeshElm2VtxConn(MPMesh_ptr p_mpmesh, int size1, int size2, int* array);
void polympo_setMeshElm2ElmConn(MPMesh_ptr p_mpmesh, int size1, int size2, int* array);

//Mesh fields
void polympo_setMeshVelArray(MPMesh_ptr p_mpmesh, int size, double* array);
void polympo_getMeshVelArray(MPMesh_ptr p_mpmesh, int size, double* array);
}
#endif
