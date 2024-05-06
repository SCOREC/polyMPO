#pragma once

#include "pmpo_utils.hpp"

#define MP_ACTIVE 1
#define MP_DELETE -1

typedef void* MPMesh_ptr;
//Function that receives void* and returns an int
typedef int (*IntVoidFunc)(void*);

using space_t = Kokkos::DefaultExecutionSpace::memory_space;

/**
 * Attention: this typedef is LayoutLeft, meaning that the first
 * index is the contiguous one. This matches the Fortran and GPU conventions for
 * allocations.
 */
//TODO: order of these typedefs to be done later
template<typename DataT>
using kkViewHostU = Kokkos::View<
          DataT,
          Kokkos::LayoutLeft,
          Kokkos::DefaultHostExecutionSpace,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

typedef kkViewHostU<double*> kkDblViewHostU;//TODO:put it to mesh.hpp             
typedef kkViewHostU<polyMPO::vec2d_t*> kkVec2dViewHostU;//TODO:put it to mesh.hpp
typedef kkViewHostU<double**> kkDbl2dViewHostU;//TODO:put it somewhere else (maybe)
typedef kkViewHostU<int**> kkInt2dViewHostU;//TODO:put it somewhere else (maybe)
typedef kkViewHostU<int*> kkIntViewHostU;//TODO:put it somewhere else (maybe)

template <typename DataT>
auto create_mirror_view_and_copy(DataT array, const int size){
  kkViewHostU<DataT> temp_host(array, size);
  return Kokkos::create_mirror_view_and_copy(space_t(), temp_host);
}
