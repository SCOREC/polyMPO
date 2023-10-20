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
void polympo_createMPs(MPMesh_ptr p_mpmesh, int numElms, int numMPs, int* mpsPerElm, int* mp2Elm, int* isMPActive);
void polympo_getMPCurElmID(MPMesh_ptr p_mpmesh, int numMPs, int* elmIDs);

//MP slices
void polympo_getMPPositions_f(MPMesh_ptr p_mpmesh, int numComps, int numMPs, double* mpPositionsIn);
void polympo_setMPPositions(MPMesh_ptr p_mpmesh, int numComps, int numMPs, double* mpPositionsIn);
void polympo_setMPVel(MPMesh_ptr p_mpmesh, int size, double* array);
void polympo_getMPVel(MPMesh_ptr p_mpmesh, int size, double* array);

//Mesh info
void polympo_startMeshFill(MPMesh_ptr p_mpmesh);
void polympo_endMeshFill(MPMesh_ptr p_mpmesh);
void polympo_checkMeshMaxSettings(MPMesh_ptr p_mpmesh, int maxEdges, int vertexDegree);
void polympo_setMeshTypeGeneralPoly(MPMesh_ptr p_mpmesh);
void polympo_setMeshTypeCVTPoly(MPMesh_ptr p_mpmesh);
void polympo_setMeshGeomTypePlanar(MPMesh_ptr p_mpmesh);
void polympo_setMeshGeomTypeSpherical(MPMesh_ptr p_mpmesh);
void polympo_setMeshSphereRadius(MPMesh_ptr p_mpmesh, double sphereRadius);
void polympo_setMeshNumVtxs(MPMesh_ptr p_mpmesh, int numVtxs);
int polympo_getMeshNumVtxs(MPMesh_ptr p_mpmesh);
void polympo_setMeshNumElms(MPMesh_ptr p_mpmesh, int numElms);
int polympo_getMeshNumElms(MPMesh_ptr p_mpmesh);
void polympo_setMeshNumEdgesPerElm(MPMesh_ptr p_mpmesh, int nCells, int* array);
void polympo_setMeshElm2VtxConn(MPMesh_ptr p_mpmesh, int maxEdges, int nCells, int* array);
void polympo_setMeshElm2ElmConn(MPMesh_ptr p_mpmesh, int maxEdges, int nCells, int* array);

//Mesh fields
void polympo_setMeshVtxCoords(MPMesh_ptr p_mpmesh, int nVertices, double* xArray, double* yArray, double* zArray);
void polympo_getMeshVtxCoords(MPMesh_ptr p_mpmesh, int nVertices, double* xArray, double* yArray, double* zArray);
void polympo_setMeshOnSurfVeloIncr(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array);//vec2d
void polympo_getMeshOnSurfVeloIncr(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array);//vec2d
void polympo_setMeshOnSurfDispIncr(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array);//vec2d
void polympo_getMeshOnSurfDispIncr(MPMesh_ptr p_mpmesh, int nComps, int nVertices, double* array);//vec2d
}
#endif
