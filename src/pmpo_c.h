#ifndef PMO_C_H
#define PMO_C_H

#include <mpi.h>
#include "pmpo_defines.h"

extern "C" {

//initialize and finalize
void polympo_initialize_f();
void polympo_finalize_f();

//create/delete MpMesh object
MPMesh_ptr polympo_createMPMesh_f(const int setMeshOption, const int setMPOption);
void polympo_deleteMPMesh_f(MPMesh_ptr p_mpmesh);

//set MPI communicator
void polympo_setMPICommunicator_f(MPMesh_ptr p_mpmesh, MPI_Fint fcomm);
//TODO: add a function to get communicator
void polympo_startCommunication_f(MPMesh_ptr p_mpmesh);

//MP info
void polympo_createMPs_f(MPMesh_ptr p_mpmesh, const int numElms, const int numMPs, int* mpsPerElm, const int* mp2Elm, const int* isMPActive);
void polympo_startRebuildMPs_f(MPMesh_ptr p_mpmesh, const int numMPs, const int* allTgtMpElmIn, const int* addedMPMask);
void polympo_startRebuildMPs2_f(MPMesh_ptr p_mpmesh, const int size1, const int* arg1, const int size2, const int size3, int* arg2, int* arg3);
int polympo_getMPCount_f(MPMesh_ptr p_mpmesh);
void polympo_finishRebuildMPs_f(MPMesh_ptr p_mpmesh);
void polympo_setAppIDFunc_f(MPMesh_ptr p_mpmesh, IntVoidFunc getNext, void* appIDs);
void polympo_getMPTgtElmID_f(MPMesh_ptr p_mpmesh, const int numMPs, int* elmIDs);
void polympo_getMPCurElmID_f(MPMesh_ptr p_mpmesh, const int numMPs, int* elmIDs);
void polympo_setMPLatLonRotatedFlag_f(MPMesh_ptr p_mpmesh, const int isRotateFlag);

//MP slices
//Positions
void polympo_setMPPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpPositionsIn);
void polympo_getMPPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpPositionsIn);
void polympo_setMPTgtPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpPositionsIn);
void polympo_getMPTgtPositions_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpPositionsIn);
//LatLon
void polympo_setMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpRotLatLonIn);
void polympo_getMPRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpRotLatLonHost);
void polympo_setMPTgtRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpRotLatLonIn);
void polympo_getMPTgtRotLatLon_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpRotLatLonHost);
//MP fields
void polympo_setMPMass_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpMassIn);
void polympo_getMPMass_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpMassHost);
void polympo_setMPVel_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpVelIn);
void polympo_getMPVel_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpVelHost);
void polympo_calculateMPStrainRate_f(MPMesh_ptr p_mpmesh);
void polympo_setMPStrainRate_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpStrainRateIn);
void polympo_getMPStrainRate_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpStrainRateHost);
void polympo_calculateMPStress_f(MPMesh_ptr p_mpmesh, const int constitutive_model);
void polympo_setMPStress_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, const double* mpStressIn);
void polympo_getMPStress_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* mpStressHost);
void polympo_setAreaMP_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* areaMPHost);
void polympo_setIcePressureMP_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* icePressureMPHost);
void polympo_getReplacementPressureMP_f(MPMesh_ptr p_mpmesh, const int nComps, const int numMPs, double* replacementPressureMPHost);

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
void polympo_setMeshNumVtxsOwned_f(MPMesh_ptr p_mpmesh, const int numVtxsOwned);
void polympo_getMeshNumVtxs_f(MPMesh_ptr p_mpmesh, int & numVtxs);
void polympo_setMeshNumElms_f(MPMesh_ptr p_mpmesh, const int numElms);
void polympo_getMeshNumElms_f(MPMesh_ptr p_mpmesh, int & numElms);
void polympo_setMeshNumEdgesPerElm_f(MPMesh_ptr p_mpmesh, const int nCells, const int* array);
void polympo_setMeshElm2VtxConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, const int* array);
void polympo_setMeshElm2ElmConn_f(MPMesh_ptr p_mpmesh, const int maxEdges, const int nCells, const int* array);
void polympo_setOwningProc_f(MPMesh_ptr p_mpmesh, const int nCells, const int* array);
void polympo_setOwningProcVertex_f(MPMesh_ptr p_mpmesh, const int nVertices, const int* array);
void polympo_setElmGlobal_f(MPMesh_ptr p_mpmesh, const int nCells, const int* array);
void polympo_setVtxGlobal_f(MPMesh_ptr p_mpmesh, const int nVertices, const int* array);
void polympo_setInteriorVertex_f(MPMesh_ptr p_mpmesh, const int nVertices, const int* array);

//Mesh fields
int polympo_getMeshFVtxType_f();
int polympo_getMeshFElmType_f();
void polympo_setMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* xArray, const double* yArray, const double* zArray);
void polympo_getMeshVtxCoords_f(MPMesh_ptr p_mpmesh, const int nVertices, double* xArray, double* yArray, double* zArray);
void polympo_setMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* latitude);
void polympo_getMeshVtxRotLat_f(MPMesh_ptr p_mpmesh, const int nVertices, double* latitude);
void polympo_setMeshVtxVel_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* uVelocity, const double* vVelocity);
void polympo_getMeshVtxVel_f(MPMesh_ptr p_mpmesh, const int nVertices, double* uVelocity, double* vVelocity);
void polympo_setMeshVtxMass_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* vtxMass);
void polympo_getMeshVtxMass_f(MPMesh_ptr p_mpmesh, const int nVertices, double* vtxMass);
void polympo_setMeshElmMass_f(MPMesh_ptr p_mpmesh, const int nCells, const double* elmMass);
void polympo_getMeshElmMass_f(MPMesh_ptr p_mpmesh, const int nCells, double* elmMass);
void polympo_setMeshVtxOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array);//vec2d
void polympo_getMeshVtxOnSurfVeloIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array);//vec2d
void polympo_setMeshVtxOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, const double* array);//vec2d
void polympo_getMeshVtxOnSurfDispIncr_f(MPMesh_ptr p_mpmesh, const int nComps, const int nVertices, double* array);//vec2d
void polympo_setMeshElmCenter_f(MPMesh_ptr p_mpmesh, const int nCells, const double* xArray, const double* yArray, const double* zArray);
void polympo_getMeshElmCenter_f(MPMesh_ptr p_mpmesh, const int nCells, double* xArray, double* yArray, double* zArray);
void polympo_setMeshDualTriangleArea_f(MPMesh_ptr p_mpmesh, const int nVertices, const double* areaTriangle);
void polympo_getMeshDualTriangleArea_f(MPMesh_ptr p_mpmesh, const int nVertices, double* areaTriangle);
void polympo_setGnomonicProjection_f(MPMesh_ptr p_mpmesh);
void polyMPO_setTanLatVertexRotatedOverRadius_f(MPMesh_ptr p_mpmesh, const int nVertices, double* array); 
void polympo_setElasticTimeStep_f(MPMesh_ptr p_mpmesh, const double elasticTimeStep);
void polympo_setDynamicTimeStep_f(MPMesh_ptr p_mpmesh, const double dynamicTimeStep);
void polympo_setSolveStressMesh_f(MPMesh_ptr p_mpmesh, const int nCells, int* array);
void polympo_setSolveVelocityMesh_f(MPMesh_ptr p_mpmesh, const int nVertices, int* array);
void polympo_calculateStressDivergence_f(MPMesh_ptr p_mpmesh);
void polympo_getStressDivergence_f(MPMesh_ptr p_mpmesh, const int nVertices, double* uArray, double* vArray);
void polympo_setTotalMassVtx_f(MPMesh_ptr p_mpmesh, const int nVertices, double* array);
void polympo_set_airStress_f(MPMesh_ptr mpmesh, const int nVertices, double* uArray,double* vArray);
void polympo_set_surfaceTiltForce_f(MPMesh_ptr p_mpmesh, const int nVertices, double* uArray, double* vArray);
void polympo_set_totalMassVertexfVertex_f(MPMesh_ptr p_mpmesh, const int nVertices, double* array);
void polympo_set_oceanStress_f(MPMesh_ptr p_mpmesh, const int nVertices, double* uArray, double* vArray);
void polympo_set_oceanStressCoefficient_f(MPMesh_ptr p_mpmesh, const int nVertices, double* array);
void polympo_velocity_grid_solve_f(MPMesh_ptr p_mpmesh);

// Advection calculations
void polympo_push_f(MPMesh_ptr p_mpmesh);
void polympo_push_ahead_f(MPMesh_ptr p_mpmesh);
bool polympo_push1P_f(MPMesh_ptr p_mpmesh);
void polympo_push_swap_f(MPMesh_ptr p_mpmesh);
void polympo_push_swap_pos_f(MPMesh_ptr p_mpmesh);

// Reconstruction of variables from MPs to mesh vertices
void polympo_setReconstructionOfMass_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType);
void polympo_setReconstructionOfVel_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType);
void polympo_setReconstructionOfStrainRate_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType);
void polympo_setReconstructionOfStress_f(MPMesh_ptr p_mpmesh, const int order, const int meshEntType);
void polympo_applyReconstruction_f(MPMesh_ptr p_mpmesh);
//Simpler/cleaner way of calling reconstruction
void polympo_reconstruct_coeff_with_MPI_f(MPMesh_ptr p_mpmesh);
void polympo_reconstruct_iceArea_with_MPI_f(MPMesh_ptr p_mpmesh);
void polympo_reconstruct_velocity_with_MPI_f(MPMesh_ptr p_mpmesh);

void polympo_init_deluDyn_f(MPMesh_ptr p_mpmesh);
void polympo_aggregate_deluDyn_f(MPMesh_ptr p_mpmesh);

// Timing
void polympo_enableTiming_f();
void polympo_summarizeTime_f();
void polympo_setTimingVerbosity_f(int v);
}
#endif
