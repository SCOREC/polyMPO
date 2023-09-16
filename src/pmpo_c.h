#ifndef PMO_C_H
#define PMO_C_H

#include <mpi.h>

extern "C" {

typedef void* MPMesh_ptr;

//initialize and finalize
void polympo_initialize_f();
void polympo_finalize_f();

//create/delete MpMesh object
MPMesh_ptr polympo_createMPMesh_f(int setMeshOption, int setMPOption);
void polympo_deleteMPMesh_f(MPMesh_ptr p_mpmesh);

//set MPI communicator
void polympo_setMPICommunicator_f(MPI_Fint fcomm);//TODO:is MPI_Fint best? or something else
//TODO: add a function to get communicator

//MP info
//TODO 

//MP slices
void polympo_setMPVelArray_f(MPMesh_ptr p_mpmesh, int size, double* array);
void polympo_getMPVelArray_f(MPMesh_ptr p_mpmesh, int size, double* array);

//Mesh info
void polympo_checkMeshMaxSettings_f(MPMesh_ptr p_mpmesh, int maxEdges, int vertexDegree);
void polympo_setMeshTypeGeneralPoly_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshTypeCVTPoly_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshGeomTypePlanar_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshGeomTypeSpherical_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshSphereRadius_f(MPMesh_ptr p_mpmesh, double sphereRadius);
void polympo_setMeshNumVtxs_f(MPMesh_ptr p_mpmesh, int numVtxs);
void polympo_setMeshNumElms_f(MPMesh_ptr p_mpmesh, int numElms);
void polympo_setMeshVtxCoords_f(MPMesh_ptr p_mpmesh, int nVertices, double* xArray, double* yArray, double* zArray);
void polympo_setMeshNumEdgesPerElm_f(MPMesh_ptr p_mpmesh, int nCells, int* array);
void polympo_setMeshElm2VtxConn_f(MPMesh_ptr p_mpmesh, int maxEdges, int nCells, int* array);
void polympo_setMeshElm2ElmConn_f(MPMesh_ptr p_mpmesh, int maxEdges, int nCells, int* array);

//Mesh fields
void polympo_setMeshVelArray_f(MPMesh_ptr p_mpmesh, int size, double* array);
void polympo_getMeshVelArray_f(MPMesh_ptr p_mpmesh, int size, double* array);
}
#endif
