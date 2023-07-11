#ifndef PMO_C_H
#define PMO_C_H

#include <mpi.h>

/**
 * Attention: these typedefs is LayoutLeft, meaning that the first
 * index is the contiguous one. This matches the Fortran and GPU conventions for
 * allocations.
 */
typedef Kokkos::View< double*[vec2d_nEntries],
                      Kokkos::LayoutLeft,
                      Kokkos::DefaultHostExecutionSpace,
                      Kokkos::MemoryTraits<Kokkos::Unmanaged>
                    > kkDblViewHostU;
typedef Kokkos::View< int**,
                      Kokkos::LayoutLeft,
                      Kokkos::DefaultHostExecutionSpace,
                      Kokkos::MemoryTraits<Kokkos::Unmanaged>
                    > kkInt2DViewHostU;

extern "C" {

typedef void* mpmesh;

void polympo_initialize();
mpmesh polympo_createMpMesh();
void polympo_deleteMpMesh(mpmesh mpMesh);
void polympo_finalize();

void polympo_setMPVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_getMPVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_setCommunicator(MPI_Fint fcomm);

void polympo_setMeshVelArray(mpmesh mpMeshIn, int size, double* array);
void polympo_getMeshVelArray(mpmesh mpMeshIn, int size, double* array);

void polympo_setMeshVtxCoords(mpmesh mpMeshIn, int size, double* xArray, double* yArray, double* zArray);
void polympo_setMeshElm2VtxConn(mpmesh mpMeshIn, int size, int** array);
void polympo_setMeshVtx2ElmConn(mpmesh mpMeshIn, int size, int** array);
void polympo_setMeshElm2ElmConn(mpmesh mpMeshIn, int size, int** array);
}
#endif
