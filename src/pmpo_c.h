#ifndef PMO_C_H
#define PMO_C_H

#include <mpi.h>

extern "C" {

typedef void* MPMesh_ptr;

//initialize and finalize
void polympo_initialize_f();
void polympo_finalize_f();

//create/delete MpMesh object
MPMesh_ptr polympo_createMPMesh_f(const int setMeshOption, const int setMPOption);
void polympo_deleteMPMesh_f(MPMesh_ptr p_mpmesh);

//set MPI communicator
void polympo_setMPICommunicator_f(MPI_Fint fcomm);//TODO:is MPI_Fint best? or something else
//TODO: add a function to get communicator

//MP info
void polympo_createMPs_f(MPMesh_ptr p_mpmesh, const int numElms, const int numMPs, int* mpsPerElm, const int* mp2Elm, const int* isMPActive);
void polympo_startRebuildMPs_f(MPMesh_ptr p_mpmesh, const int numMPs, const int* allTgtMpElmIn, const int* addedMPMask);
void polympo_finishRebuildMPs_f(MPMesh_ptr p_mpmesh);
void polympo_getMPCurElmID_f(MPMesh_ptr p_mpmesh, const int numMPs, int* elmIDs);
void polympo_setMPLatLonRotatedFlag_f(MPMesh_ptr p_mpmesh, const int isRotateFlag);

//MP slices
void polympo_setMPPositions_f(MPMesh_ptr p_mpmesh, const int numComps, const int numMPs, const double* mpPositionsIn);
void polympo_getMPPositions_f(MPMesh_ptr p_mpmesh, const int numComps, const int numMPs, double* mpPositionsIn);
void polympo_setMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int numComps, const int numMPs, const double* mpRotLatLonIn);
void polympo_getMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int numComps, const int numMPs, double* mpRotLatLonHost);
void polympo_setMPVel_f(MPMesh_ptr p_mpmesh, const int size, const double* array);
void polympo_getMPVel_f(MPMesh_ptr p_mpmesh, const int size, double* array);

//Mesh info
void polympo_startMeshFill_f(MPMesh_ptr p_mpmesh);
void polympo_endMeshFill_f(MPMesh_ptr p_mpmesh);
void polympo_checkMeshMaxSettings_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int vertexDegree);
void polympo_setMeshTypeGeneralPoly_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshTypeCVTPoly_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshGeomTypePlanar_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshGeomTypeSpherical_f(MPMesh_ptr p_mpmesh);
void polympo_setMeshSphereRadius_f(MPMesh_ptr p_mpmesh, const double sphereRadius);
void polympo_setMeshNumVtxs_f(MPMesh_ptr p_mpmesh, const int numVtxs);
void polympo_getMeshNumVtxs_f(MPMesh_ptr p_mpmesh, int & numVtxs);
void polympo_setMeshNumElms_f(MPMesh_ptr p_mpmesh, const int numElms);
void polympo_getMeshNumElms_f(MPMesh_ptr p_mpmesh, int & numElms);
void polympo_setMeshNumEdgesPerElm_f(MPMesh_ptr p_mpmesh, const int nCells, const int* array);
void polympo_setMeshElm2VtxConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, const int* array);
void polympo_setMeshElm2ElmConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, const int* array);

//Mesh fields
void polympo_setMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* xArray, const double* yArray, const double* zArray);
void polympo_getMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, double* xArray, double* yArray, double* zArray);
void polympo_setMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* latitude);
void polympo_getMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, double* latitude);
void polympo_setMeshOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array);//vec2d
void polympo_getMeshOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array);//vec2d
void polympo_setMeshOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array);//vec2d
void polympo_getMeshOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array);//vec2d

// calculations
void polympo_push_f(MPMesh_ptr p_mpmesh);
}
#endif
