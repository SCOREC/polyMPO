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

template <typename DataT>
auto create_mirror_view_and_copy(DataT array, const int size){
  kkViewHostU<DataT> temp_host(array, size);
  return Kokkos::create_mirror_view_and_copy(space_t(), temp_host);
}
