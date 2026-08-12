#ifndef POLYMPO_MPM_H
#define POLYMPO_MPM_H

#include "pmpo_utils.hpp"
#include "pmpo_mesh.hpp"
#include "pmpo_materialPoints.hpp"
#include <cstdlib>
#include <iostream>
#include <utility>
#include <stdexcept>
#include <string>
#ifdef KOKKOS_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

// ----------------------------------------------------------------------
// Halo-exchange optimization: pure derived-datatype send (pack path fully
// removed) + MPI_Waitany incremental unpack.
//
// communicate_and_take_halo_contributions1_improved() (CUDA_AWARE_MPI path
// only - the CPU-staged fallback below is untouched):
//
//   1. Eliminate the pack step / sendDataGPU buffer, unconditionally. Each
//      neighbor proc gets its own committed MPI derived datatype
//      (MPI_Type_create_hindexed_block over a "rowType" built from
//      meshField's own strides) that gathers that proc's rows directly out
//      of meshField's GPU memory. MPI_Isend reads straight from the field -
//      no pack kernel, no extra buffer, no fallback to a pack path.
//
//   2. Hide the unpack step behind MPI_Waitall latency instead of paying
//      for it serially afterward. Recv requests are drained with
//      MPI_Waitany instead of MPI_Waitall, and each neighbor's contribution
//      is unpacked (async kernel launch, no intermediate fence) the instant
//      that neighbor's message lands, so unpacking early arrivals overlaps
//      with waiting on stragglers. A single Kokkos::fence() at the end
//      guarantees every launched unpack kernel has completed before
//      meshField is used downstream.
//
// The receive side still uses a batched contiguous GPU buffer + a real
// unpack kernel (not a derived receive-datatype): for op==0 the unpack is a
// scatter-ADD (Kokkos::atomic_add), because a single owned vertex can
// receive contributions from multiple different remote ranks across
// separate messages, and a derived datatype can only place bytes at an
// offset, not reduce concurrently-arriving values. Only the gather (send)
// side is a pure "pick these rows" operation, which is what derived
// datatypes are actually good for here.
//
// Known risk, confirmed relevant on this codebase's field layout
// (LayoutLeft): sending directly from meshField's memory means MPI's
// datatype engine has to gather scattered elements (MPI_Type_vector nested
// in MPI_Type_create_hindexed_block) instead of a few large contiguous
// blocks, which is expensive for many MPI implementations' datatype
// engines to execute well relative to a hand-written parallel GPU pack
// kernel. It also means MPI is handed a pointer straight into meshField's
// own (Kokkos-managed, possibly CudaMallocAsync pool-allocated) memory,
// rather than a dedicated raw-cudaMalloc'd staging buffer - see
// RawCudaMPIBuffer's comment below for why that distinction matters on
// Cray MPICH/GTL. A pre-change, fully protected copy of this file's
// send-side logic is preserved in pmpo_MPMesh_backup.hpp and in git history
// (the "Add CUDA-aware MPI halo exchange" / "Add CUDA-aware communication
// path to MPMesh" commits) for comparison/revert.
// ----------------------------------------------------------------------

namespace polyMPO{

template <MeshFieldIndex>
const MaterialPointSlice meshFieldIndexToMPSlice;
template <> const MaterialPointSlice meshFieldIndexToMPSlice < MeshF_Vel            > = MPF_Vel;
template <> const MaterialPointSlice meshFieldIndexToMPSlice < MeshF_VtxMass        > = MPF_Mass;
template <> const MaterialPointSlice meshFieldIndexToMPSlice < MeshF_ElmMass        > = MPF_Mass;
template <> const MaterialPointSlice meshFieldIndexToMPSlice < MeshF_RotLatLonIncr  > = MPF_Rot_Lat_Lon_Incr;
template <> const MaterialPointSlice meshFieldIndexToMPSlice < MeshF_OnSurfVeloIncr > = MPF_Vel_Incr;

#define maxMPsPerElm 8

class MPMesh{

  public:

    Mesh* p_mesh;
    MaterialPoints* p_MPs;

    //For MPI Communication
    int numOwnersTot, numHalosTot;
    std::vector<int> numOwnersOnOtherProcs;
    std::vector<int> numHalosOnOtherProcs;
    std::vector<int> haloOwnerProcs;
    std::vector<std::vector<int>> haloOwnerLocalIDs;
    std::vector<std::vector<int>> ownerOwnerLocalIDs;
    std::vector<std::vector<int>> ownerHaloLocalIDs;

    void startCommunication();

    void communicate_and_take_halo_contributions(
        const Kokkos::View<double**>& meshField,
        int nEntities,
        int numEntries,
        int mode,
        int op);

    // Original CPU-staging function
    template <typename ViewType>
    void communicate_and_take_halo_contributions1(
        const ViewType& meshField,
        int nEntities,
        int numEntries,
        int mode,
        int op){

      int self;
      MPI_Comm comm = p_MPs->getMPIComm();
      MPI_Comm_rank(comm, &self);

      Kokkos::Timer timer;
      auto reconVals_host =
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), meshField);

      pumipic::RecordTime("SD: GPU-CPU copy-" + std::to_string(self), timer.seconds());

      timer.reset();
      std::vector<std::vector<int>> recvIDVec;
      std::vector<std::vector<double>> recvDataVec;

      pumipic::RecordTime("SD: Recv Vec Allocation-" + std::to_string(self), timer.seconds());

      timer.reset();

      communicateFields1(
          reconVals_host,
          nEntities,
          numEntries,
          mode,
          recvIDVec,
          recvDataVec);

      pumipic::RecordTime("SD: IP Comm-" + std::to_string(self), timer.seconds());

      timer.reset();

      int numProcsTot = recvIDVec.size();

      int totalSize = 0;
      std::vector<int> offsets(numProcsTot, 0);

      for(int i = 0; i < numProcsTot; i++){
        offsets[i] = totalSize;
        totalSize += recvIDVec[i].size();
      }

      Kokkos::View<int*> recvIDGPU("recvIDGPU", totalSize);
      auto hostView =
          Kokkos::View<int*, Kokkos::HostSpace>("recvIDCPU", totalSize);

      for(int i = 0; i < numProcsTot; i++){
        std::copy(
            recvIDVec[i].begin(),
            recvIDVec[i].end(),
            hostView.data() + offsets[i]);
      }

      pumipic::RecordTime("SD: Flatten IDs-" + std::to_string(self), timer.seconds());

      timer.reset();

      Kokkos::deep_copy(recvIDGPU, hostView);
      Kokkos::fence();

      pumipic::RecordTime("SD: Copy CPU-GPU-" + std::to_string(self), timer.seconds());

      timer.reset();

      int totalSize_data = 0;
      std::vector<int> offsets_data(numProcsTot, 0);

      for(int i = 0; i < numProcsTot; i++){
        offsets_data[i] = totalSize_data;
        totalSize_data += recvDataVec[i].size();
      }

      Kokkos::View<double*> recvDataGPU("recvDataGPU", totalSize_data);
      auto hostView_data =
          Kokkos::View<double*, Kokkos::HostSpace>("recvDataCPU", totalSize_data);

      for(int i = 0; i < numProcsTot; i++){
        std::copy(
            recvDataVec[i].begin(),
            recvDataVec[i].end(),
            hostView_data.data() + offsets_data[i]);
      }

      pumipic::RecordTime("SD: Flatten Data-" + std::to_string(self), timer.seconds());

      timer.reset();

      Kokkos::deep_copy(recvDataGPU, hostView_data);
      Kokkos::fence();

      assert(totalSize_data == totalSize * numEntries);

      for(int i = 0; i < numProcsTot; i++){
        assert(recvDataVec[i].size() == recvIDVec[i].size() * numEntries);
      }

      pumipic::RecordTime("SD: Copy CPU-GPU2-" + std::to_string(self), timer.seconds());

      timer.reset();

      if(op == 0){
        Kokkos::parallel_for(
            "halo contribution add",
            recvIDGPU.size(),
            KOKKOS_LAMBDA(const int i){
              const int vertex = recvIDGPU(i);

              for(int k = 0; k < numEntries; k++){
                Kokkos::atomic_add(
                    &meshField(vertex,k),
                    recvDataGPU(i * numEntries + k));
              }
            });
      }
      else{
        Kokkos::parallel_for(
            "halo contribution assign",
            recvIDGPU.size(),
            KOKKOS_LAMBDA(const int i){
              const int vertex = recvIDGPU(i);

              for(int k = 0; k < numEntries; k++){
                meshField(vertex,k) =
                    recvDataGPU(i * numEntries + k);
              }
            });
      }

      Kokkos::fence();

      pumipic::RecordTime("SD: Contribution" + std::to_string(self), timer.seconds());
    }


    void communicateFields(
        const std::vector<std::vector<double>>& fieldData,
        const int numEntities,
        const int numEntries,
        int mode,
        std::vector<std::vector<int>>& recvIDVec,
        std::vector<std::vector<double>>& recvDataVec);


    template <class ViewType>
    void communicateFields1(
        const ViewType& fieldData,
        const int numEntities,
        const int numEntries,
        int mode,
        std::vector<std::vector<int>>& recvIDVec,
        std::vector<std::vector<double>>& recvDataVec){

      int self, numProcsTot;

      MPI_Comm comm = p_MPs->getMPIComm();

      MPI_Comm_rank(comm, &self);
      MPI_Comm_size(comm, &numProcsTot);

      assert(numEntities == numOwnersTot + numHalosTot);

      std::vector<std::vector<double>> sendDataVec(numProcsTot);

      recvIDVec.resize(numProcsTot);
      recvDataVec.resize(numProcsTot);

      for(int i = 0; i < numProcsTot; i++){
        if(i == self) continue;

        int numToSend = 0;
        int numToRecv = 0;

        if(mode == 0){
          numToSend = numOwnersOnOtherProcs[i];
          numToRecv = numHalosOnOtherProcs[i];
        }
        else{
          numToSend = numHalosOnOtherProcs[i];
          numToRecv = numOwnersOnOtherProcs[i];
        }

        if(numToSend > 0){
          sendDataVec[i].reserve(numToSend * numEntries);
        }

        if(numToRecv > 0){
          recvDataVec[i].resize(numToRecv * numEntries);
          recvIDVec[i].resize(numToRecv);
        }
      }

      if(mode == 0){
        for(int iEnt = 0; iEnt < numHalosTot; iEnt++){
          auto ownerProc = haloOwnerProcs[iEnt];

          for(int iDouble = 0; iDouble < numEntries; iDouble++){
            sendDataVec[ownerProc].push_back(
                fieldData(numOwnersTot + iEnt, iDouble));
          }
        }
      }
      else if(mode == 1){
        for(size_t iProc = 0; iProc < ownerOwnerLocalIDs.size(); iProc++){
          for(auto& ownerID : ownerOwnerLocalIDs[iProc]){
            for(int iDouble = 0; iDouble < numEntries; iDouble++){
              sendDataVec[iProc].push_back(
                  fieldData(ownerID, iDouble));
            }
          }
        }
      }

      std::vector<MPI_Request> requests;
      requests.reserve(4 * numProcsTot);

      for(int proc = 0; proc < numProcsTot; proc++){
        if(proc == self) continue;

        if(mode == 0 && numHalosOnOtherProcs[proc]){
          assert(recvIDVec[proc].size() ==
                 static_cast<size_t>(numHalosOnOtherProcs[proc]));

          assert(recvDataVec[proc].size() ==
                 recvIDVec[proc].size() * static_cast<size_t>(numEntries));

          MPI_Request req3;
          MPI_Request req4;

          MPI_Irecv(
              recvIDVec[proc].data(),
              recvIDVec[proc].size(),
              MPI_INT,
              proc,
              1,
              comm,
              &req3);

          MPI_Irecv(
              recvDataVec[proc].data(),
              recvDataVec[proc].size(),
              MPI_DOUBLE,
              proc,
              2,
              comm,
              &req4);

          requests.push_back(req3);
          requests.push_back(req4);
        }

        if(mode == 0 && numOwnersOnOtherProcs[proc]){
          assert(haloOwnerLocalIDs[proc].size() ==
                 static_cast<size_t>(numOwnersOnOtherProcs[proc]));

          assert(sendDataVec[proc].size() ==
                 haloOwnerLocalIDs[proc].size() * static_cast<size_t>(numEntries));

          MPI_Request req1;
          MPI_Request req2;

          MPI_Isend(
              haloOwnerLocalIDs[proc].data(),
              haloOwnerLocalIDs[proc].size(),
              MPI_INT,
              proc,
              1,
              comm,
              &req1);

          MPI_Isend(
              sendDataVec[proc].data(),
              sendDataVec[proc].size(),
              MPI_DOUBLE,
              proc,
              2,
              comm,
              &req2);

          requests.push_back(req1);
          requests.push_back(req2);
        }

        if(mode == 1 && numOwnersOnOtherProcs[proc]){
          MPI_Request req3;
          MPI_Request req4;

          MPI_Irecv(
              recvIDVec[proc].data(),
              recvIDVec[proc].size(),
              MPI_INT,
              proc,
              1,
              comm,
              &req3);

          MPI_Irecv(
              recvDataVec[proc].data(),
              recvDataVec[proc].size(),
              MPI_DOUBLE,
              proc,
              2,
              comm,
              &req4);

          requests.push_back(req3);
          requests.push_back(req4);
        }

        if(mode == 1 && numHalosOnOtherProcs[proc]){
          MPI_Request req1;
          MPI_Request req2;

          MPI_Isend(
              ownerHaloLocalIDs[proc].data(),
              ownerHaloLocalIDs[proc].size(),
              MPI_INT,
              proc,
              1,
              comm,
              &req1);

          MPI_Isend(
              sendDataVec[proc].data(),
              sendDataVec[proc].size(),
              MPI_DOUBLE,
              proc,
              2,
              comm,
              &req2);

          requests.push_back(req1);
          requests.push_back(req2);
        }
      }

      MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    }


    MPMesh(Mesh* inMesh, MaterialPoints* inMPs):
      p_mesh(inMesh),
      p_MPs(inMPs) {
    };


    ~MPMesh(){
      delete p_mesh;
      delete p_MPs;
    }


    // MP advection and tracking
    void CVTTrackingEdgeCenterBased(Vec2dView dx);
    void CVTTrackingElmCenterBased(const int printVTPIndex = -1);
    void T2LTracking(Vec2dView dx);
    bool push1P();
    void push_ahead();
    void push_swap();
    void push_swap_pos();
    void push();


    // Used before advection to interpolate fields from mesh to MPs
    // And also before reconstruction
    void calcBasis();


    // Reconstruction
    DoubleView assemblyV0();

    template <MeshFieldIndex meshFieldIndex>
    void assemblyVtx0();

    template <MeshFieldIndex meshFieldIndex>
    void assemblyElm0();

    template <MeshFieldIndex meshFieldIndex>
    void assemblyVtx1();

    void reconstruct_coeff_full();

    void invertMatrix(
        const Kokkos::View<double**>& vtxMatrices,
        const double& radius);

    Kokkos::View<double*[vec3d_nEntries][vec4d_nEntries]> precomputedVtxCoeffs_new;
    Kokkos::View<double*> nearAnEdge;
    Kokkos::View<double*> vtxMatrixMass;


    // Not used currently
    std::map<MeshFieldIndex, std::function<void()>> reconstructSlice =
        std::map<MeshFieldIndex, std::function<void()>>();

    template <MaterialPointSlice index>
    DoubleView wtScaAssembly();

    template <MaterialPointSlice index>
    Vec2dView wtVec2Assembly();

    template <MeshFieldIndex meshFieldIndex>
    void assembly(
        int order,
        MeshFieldType type,
        bool basisWeightFlag,
        bool massWeightFlag);

    template<MeshFieldIndex meshFieldIndex>
    void setReconstructSlice(
        int order,
        MeshFieldType type);

    void reconstructSlices();

    void printVTP_mesh(int printVTPIndex);

    void writeMPTrackingVTP(
        int printVTPIndex,
        int numMPs,
        const Vec3dView& history,
        const Vec3dView& resultLeft,
        const Vec3dView& resultRight,
        const Vec3dView& mpTgtPosArray);

    void calculateStrain();
    void calculateStress(const int constitutive_relation);
    void calculateStressDivergence();



#ifdef CUDA_AWARE_MPI

    // Use explicit CUDA memory space for MPI device buffers when CUDA is
    // enabled. This avoids ambiguity in the default Kokkos::View memory space
    // and gives Cray MPICH/GTL plain CUDA allocations to register/export.
#ifdef KOKKOS_ENABLE_CUDA
    // Plain cudaMalloc'd device buffer, exposed as an unmanaged Kokkos::View
    // via .view(). Used only for the recv-side buffers below (recvIDGPU /
    // recvDataGPU) that get handed directly to MPI_Irecv under CUDA-aware
    // MPI. The send side no longer needs a plain buffer here - it sends
    // directly out of meshField via a derived datatype.
    //
    // Why not just a Kokkos::View<T*, Kokkos::CudaSpace>: when Kokkos is
    // built with Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=ON (the default since
    // Kokkos 4.2), View allocations use cudaMallocAsync/memory pools.
    // cuIpcGetMemHandle (which Cray MPICH/GTL uses for intra-node GPU-to-GPU
    // sends) rejects pool allocations with CUDA_ERROR_INVALID_VALUE. Since
    // polyMPO isn't allowed to touch the Kokkos build config, these buffers
    // bypass Kokkos's allocator entirely via a direct cudaMalloc, which
    // cuIpcGetMemHandle always accepts, regardless of how the rest of Kokkos
    // (or the rest of the app's Views) is configured.
    template <typename T>
    struct RawCudaMPIBuffer{
      T* ptr = nullptr;
      size_t count = 0;

      void allocate(size_t n){
        free();
        count = n;
        if(n > 0){
          cudaError_t err = cudaMalloc(&ptr, n * sizeof(T));
          if(err != cudaSuccess){
            throw std::runtime_error(
                std::string("RawCudaMPIBuffer: cudaMalloc failed: ") +
                cudaGetErrorString(err));
          }
        }
      }

      void free(){
        if(ptr != nullptr){ cudaFree(ptr); ptr = nullptr; }
        count = 0;
      }

      T* data() const{ return ptr; }
      size_t size() const{ return count; }

      Kokkos::View<T*, Kokkos::CudaSpace, Kokkos::MemoryUnmanaged> view() const{
        return Kokkos::View<T*, Kokkos::CudaSpace, Kokkos::MemoryUnmanaged>(
            ptr, count);
      }

      RawCudaMPIBuffer() = default;
      ~RawCudaMPIBuffer(){ free(); }

      RawCudaMPIBuffer(const RawCudaMPIBuffer&) = delete;
      RawCudaMPIBuffer& operator=(const RawCudaMPIBuffer&) = delete;

      RawCudaMPIBuffer(RawCudaMPIBuffer&& other) noexcept{
        ptr = other.ptr; count = other.count;
        other.ptr = nullptr; other.count = 0;
      }

      RawCudaMPIBuffer& operator=(RawCudaMPIBuffer&& other) noexcept{
        if(this != &other){
          free();
          ptr = other.ptr; count = other.count;
          other.ptr = nullptr; other.count = 0;
        }
        return *this;
      }
    };

    using CudaAwareMPIIntBuffer = RawCudaMPIBuffer<int>;
    using CudaAwareMPIDoubleBuffer = RawCudaMPIBuffer<double>;
#else
    // No CUDA backend: no CUDA IPC/pool-allocation concern, so just wrap a
    // normal Kokkos::View with the same .allocate()/.data()/.view()
    // interface as RawCudaMPIBuffer above, so the cache struct and its call
    // sites below don't need to branch on KOKKOS_ENABLE_CUDA.
    template <typename T>
    struct KokkosMPIBuffer{
      Kokkos::View<T*> v;

      void allocate(size_t n){
        v = Kokkos::View<T*>("cudaAwareMPIBuffer_batched", n);
      }

      T* data() const{ return v.data(); }
      size_t size() const{ return v.extent(0); }
      Kokkos::View<T*> view() const{ return v; }
    };

    using CudaAwareMPIIntBuffer = KokkosMPIBuffer<int>;
    using CudaAwareMPIDoubleBuffer = KokkosMPIBuffer<double>;
#endif

    // Cached CUDA-aware MPI communication metadata and batched GPU buffers.
    bool cudaAwareMPICacheValid = true;
    bool cudaAwareMPIDisabled = false;
    bool cudaAwareMPIEnvChecked = false;
    bool cudaAwareMPIForceCPU = false;
    bool cudaAwareMPILogged = false;

    struct CudaAwareMPIFieldCache{
      bool valid = false;
      int cachedNumProcs = -1;

      std::vector<int> sendCounts;
      std::vector<int> recvCounts;
      std::vector<int> sendOffsets; // prefix sum of sendCounts, in entities
      std::vector<int> recvOffsets; // prefix sum of recvCounts, in entities

      int totalSendCount = 0;
      int totalRecvCount = 0;

      // Receive side keeps a single batched GPU buffer (one allocation each,
      // instead of one Kokkos::View per neighbor proc). Per-proc slices are
      // [offset, offset + count) for recvIDGPU, and
      // [offset * numEntries, (offset + count) * numEntries) for
      // recvDataGPU. Still needed because op==0 unpacking is a scatter-ADD
      // (Kokkos::atomic_add, since a single owned vertex can receive
      // contributions from several different remote ranks across separate
      // messages), which a derived receive-datatype can't do on its own.
      CudaAwareMPIIntBuffer recvIDGPU;
      CudaAwareMPIDoubleBuffer recvDataGPU;

      // ---- Send side ----
      // No pack buffer/kernel anymore, unconditionally. rowType describes
      // one entity's numEntries doubles exactly as they sit in meshField's
      // own memory (built from meshField's actual strides, so it's correct
      // whether the view is LayoutRight, LayoutLeft, or padded).
      // sendProcTypes[proc] wraps rowType with that proc's list of entity
      // byte-displacements via MPI_Type_create_hindexed_block, so MPI_Isend
      // can gather directly out of meshField's GPU memory with no separate
      // pack kernel/buffer, always.
      MPI_Datatype rowType = MPI_DATATYPE_NULL;
      std::vector<MPI_Datatype> sendProcTypes; // size numProcsTot; MPI_DATATYPE_NULL where unused

      void freeTypes(){
        if(rowType != MPI_DATATYPE_NULL){
          MPI_Type_free(&rowType);
          rowType = MPI_DATATYPE_NULL;
        }

        for(auto& t : sendProcTypes){
          if(t != MPI_DATATYPE_NULL){
            MPI_Type_free(&t);
            t = MPI_DATATYPE_NULL;
          }
        }

        sendProcTypes.clear();
      }

      // Assumes this cache (and therefore the MPMesh that owns it) is torn
      // down before MPI_Finalize - same assumption the rest of this class
      // already makes about the MPI resources it holds.
      ~CudaAwareMPIFieldCache(){
        freeTypes();
      }
    };

    std::map<std::pair<int, int>, CudaAwareMPIFieldCache> cudaAwareMPICaches;

    bool cudaAwareMPIForceDisabled(){
      if(!cudaAwareMPIEnvChecked){
        const char* value = std::getenv("POLYMPO_DISABLE_CUDA_AWARE_MPI");
        cudaAwareMPIForceCPU =
            value != nullptr && value[0] != '\0' && value[0] != '0';
        cudaAwareMPIEnvChecked = true;
      }

      return cudaAwareMPIForceCPU;
    }

    // Fully CUDA-aware MPI version, batched buffer + derived-datatype send
    // variant:
    //
    // Field data is sent using GPU pointers gathered directly by a per-
    // neighbor MPI derived datatype (no pack kernel/buffer, unconditionally),
    // and received into a single batched GPU buffer per cache entry, same
    // as before. Receive IDs are cached once from the fixed halo/owner
    // mapping and are not sent every call.
    //
    // Receive completion is drained with MPI_Waitany instead of
    // MPI_Waitall: each neighbor's contribution is unpacked (async kernel
    // launch) the moment that neighbor's message arrives, so unpacking
    // already-arrived neighbors overlaps with waiting on the remaining
    // (slower) ones, instead of the previous "wait for everyone, then
    // unpack everyone" ordering. A single Kokkos::fence() at the end
    // guarantees all launched unpack kernels have completed.
    //
    // Note: an earlier version of this cache used one Kokkos::View
    // allocation per neighbor specifically to avoid handing MPI a
    // "base pointer + offset" GPU address, out of concern for CUDA IPC
    // issues on Cray MPICH/GTL. That failure mode was root-caused to
    // Kokkos allocating device Views via cudaMallocAsync (invalid for
    // cuIpcGetMemHandle), not to offset pointers themselves, and is fixed
    // by building Kokkos with -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF.
    // The same offset-pointer reasoning applies to the derived-datatype
    // send below (MPI is given meshField's base pointer plus per-entity
    // byte displacements); if you hit IPC trouble that tracks back to that,
    // the previous pack-buffer send path can be restored from
    // pmpo_MPMesh_backup.hpp / version control.
    //
    // Important:
    // This function caches communication metadata, GPU buffers, and MPI
    // derived datatypes per (mode, numEntries). If the communication
    // pattern changes, clear cudaAwareMPICaches before the next call.
    template <typename ViewType>
    void communicate_and_take_halo_contributions1_improved(
        const ViewType& meshField,
        int nEntities,
        int numEntries,
        int mode,
        int op,
        const std::string& label){

      int self, numProcsTot;

      MPI_Comm comm = p_MPs->getMPIComm();

      MPI_Comm_rank(comm, &self);
      MPI_Comm_size(comm, &numProcsTot);

      const char* diagnosticsEnv =

      std::getenv("POLYMPO_MPI_DIAGNOSTICS");

      const bool mpiDiagnostics =
      diagnosticsEnv != nullptr &&
      std::atoi(diagnosticsEnv) != 0;

      assert(mode == 0 || mode == 1);
      assert(op == 0 || op == 1);
      assert(nEntities == numOwnersTot + numHalosTot);

      if(cudaAwareMPIDisabled || cudaAwareMPIForceDisabled()){
        communicate_and_take_halo_contributions1(
            meshField,
            nEntities,
            numEntries,
            mode,
            op);
        return;
      }

#ifdef POLYMPO_VERBOSE_MPI
      if(self == 0 && !cudaAwareMPILogged){
        std::cout
            << "[CUDA_AWARE_MPI] Using batched single-buffer GPU-aware MPI path (derived-datatype send + MPI_Waitany incremental unpack) in communicate_and_take_halo_contributions1_improved()"
            << "\n";
        cudaAwareMPILogged = true;
      }
#endif

      Kokkos::Timer timer;

      if(!cudaAwareMPICacheValid){
        cudaAwareMPICaches.clear();
        cudaAwareMPICacheValid = true;
      }

      auto& cudaAwareCache =
          cudaAwareMPICaches[std::make_pair(mode, numEntries)];

      const bool needRebuild =
          (!cudaAwareCache.valid) ||
          (cudaAwareCache.cachedNumProcs != numProcsTot);

      if(needRebuild){

        // Drop any datatypes committed for a previous topology before
        // rebuilding (no-op the first time, when nothing has been built
        // yet).
        cudaAwareCache.freeTypes();

        cudaAwareCache.cachedNumProcs = numProcsTot;

        cudaAwareCache.sendCounts.assign(numProcsTot, 0);
        cudaAwareCache.recvCounts.assign(numProcsTot, 0);
        cudaAwareCache.sendOffsets.assign(numProcsTot, 0);
        cudaAwareCache.recvOffsets.assign(numProcsTot, 0);

        for(int proc = 0; proc < numProcsTot; proc++){
          if(proc == self) continue;

          if(mode == 0){
            cudaAwareCache.sendCounts[proc] = numOwnersOnOtherProcs[proc];
            cudaAwareCache.recvCounts[proc] = numHalosOnOtherProcs[proc];
          }
          else{
            cudaAwareCache.sendCounts[proc] = numHalosOnOtherProcs[proc];
            cudaAwareCache.recvCounts[proc] = numOwnersOnOtherProcs[proc];
          }
        }

        int totalSend = 0;
        int totalRecv = 0;

        for(int proc = 0; proc < numProcsTot; proc++){
          cudaAwareCache.sendOffsets[proc] = totalSend;
          totalSend += cudaAwareCache.sendCounts[proc];

          cudaAwareCache.recvOffsets[proc] = totalRecv;
          totalRecv += cudaAwareCache.recvCounts[proc];
        }

        cudaAwareCache.totalSendCount = totalSend;
        cudaAwareCache.totalRecvCount = totalRecv;

        cudaAwareCache.recvIDGPU.allocate(totalRecv);
        cudaAwareCache.recvDataGPU.allocate(totalRecv * numEntries);

        // ---- Build the send-side derived datatypes, unconditionally ----
        // Instead of a flattened host entity list + GPU pack buffer, build
        // one MPI_Datatype per neighbor proc that gathers that proc's rows
        // directly out of meshField's own memory.
        if(totalSend > 0){

          // rowType: one entity's numEntries doubles, as they actually sit
          // in meshField's memory. Query the real strides rather than
          // assuming a layout, so this is correct for LayoutRight (entries
          // within an entity contiguous - typical host default) and
          // LayoutLeft (entries within an entity strided by the entity
          // count - typical Kokkos CUDA default) alike. Note: for
          // LayoutLeft this produces an MPI_Type_vector describing
          // widely-scattered elements, which is the expensive case flagged
          // in the file-level comment above.
          const size_t strideEntity = meshField.stride_0();
          const size_t strideEntry  = meshField.stride_1();

          if(strideEntry == 1){
            MPI_Type_contiguous(numEntries, MPI_DOUBLE, &cudaAwareCache.rowType);
          }
          else{
            MPI_Type_vector(
                numEntries,
                1,
                static_cast<int>(strideEntry),
                MPI_DOUBLE,
                &cudaAwareCache.rowType);
          }
          MPI_Type_commit(&cudaAwareCache.rowType);

          std::vector<std::vector<int>> sendEntityIDsByProc(numProcsTot);

          if(mode == 0){
            for(int iEnt = 0; iEnt < numHalosTot; iEnt++){
              int ownerProc = haloOwnerProcs[iEnt];
              if(ownerProc == self) continue;

              sendEntityIDsByProc[ownerProc].push_back(numOwnersTot + iEnt);
            }
          }
          else{
            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;

              for(auto& ownerID : ownerOwnerLocalIDs[proc]){
                sendEntityIDsByProc[proc].push_back(ownerID);
              }
            }
          }

          cudaAwareCache.sendProcTypes.assign(numProcsTot, MPI_DATATYPE_NULL);

          for(int proc = 0; proc < numProcsTot; proc++){
            if(proc == self) continue;

            const int sendCount = cudaAwareCache.sendCounts[proc];
            if(sendCount <= 0) continue;

            assert(sendEntityIDsByProc[proc].size() ==
                   static_cast<size_t>(sendCount));

            std::vector<MPI_Aint> displacements(sendCount);

            for(int i = 0; i < sendCount; i++){
              displacements[i] =
                  static_cast<MPI_Aint>(sendEntityIDsByProc[proc][i]) *
                  static_cast<MPI_Aint>(strideEntity) *
                  static_cast<MPI_Aint>(sizeof(double));
            }

            MPI_Type_create_hindexed_block(
                sendCount,
                1,
                displacements.data(),
                cudaAwareCache.rowType,
                &cudaAwareCache.sendProcTypes[proc]);

            MPI_Type_commit(&cudaAwareCache.sendProcTypes[proc]);
          }
        }

        // ---- Build the flattened recv-ID list (host, then one deep_copy) ----
        if(totalRecv > 0){
          auto recvIDCPU =
              Kokkos::View<int*, Kokkos::HostSpace>(
                  "recvIDCPU_batched", totalRecv);

          if(mode == 0){
            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;

              int recvCount = cudaAwareCache.recvCounts[proc];
              if(recvCount <= 0) continue;

              assert(ownerOwnerLocalIDs[proc].size() ==
                     static_cast<size_t>(recvCount));

              int base = cudaAwareCache.recvOffsets[proc];

              for(int i = 0; i < recvCount; i++){
                recvIDCPU(base + i) = ownerOwnerLocalIDs[proc][i];
              }
            }
          }
          else{
            std::vector<int> cursor(cudaAwareCache.recvOffsets);

            for(int iEnt = 0; iEnt < numHalosTot; iEnt++){
              int ownerProc = haloOwnerProcs[iEnt];
              if(ownerProc == self) continue;

              recvIDCPU(cursor[ownerProc]) = numOwnersTot + iEnt;
              cursor[ownerProc]++;
            }

            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;

              assert(cursor[proc] ==
                     cudaAwareCache.recvOffsets[proc] +
                         cudaAwareCache.recvCounts[proc]);
            }
          }

          Kokkos::deep_copy(cudaAwareCache.recvIDGPU.view(), recvIDCPU);
        }

        cudaAwareCache.valid = true;

        if(mpiDiagnostics){
            pumipic::RecordTime(label + "_MPI_Diagnostics_CacheBuild_m" + std::to_string(mode) + "_e" + std::to_string(numEntries) +
            "_rank" + std::to_string(self), timer.seconds());
      }

        timer.reset();
      }

      // ---- No explicit pack step, unconditionally ----
      // MPI_Isend below reads meshField directly through the derived
      // datatypes built above, so there is no separate pack kernel or
      // sendDataGPU buffer to launch/fence on here. We still fence once so
      // that any of the caller's kernels which wrote meshField complete
      // before MPI starts reading it directly (mirrors the fence the old
      // pack step used to provide, just guarding meshField itself now
      // instead of a pack buffer).
      Kokkos::fence();

      if(mpiDiagnostics){
         pumipic::RecordTime(
         label + "_MPI_Diagnostics_Pack_m" + std::to_string(mode) + "_e" + std::to_string(numEntries) + "_rank" + std::to_string(self), 0.0);
    }

      timer.reset();

      double postTime = 0.0;
      double waitTime = 0.0;

      std::vector<MPI_Request> recvRequests;
      std::vector<int> recvReqProc;
      std::vector<MPI_Request> sendRequests;

      recvRequests.reserve(numProcsTot);
      recvReqProc.reserve(numProcsTot);
      sendRequests.reserve(numProcsTot);

      int mpiError = MPI_SUCCESS;

      // Data volume exchanged by this rank in this call (recorded once per
      // call, independent of the post/wait timing below), tagged by the
      // caller (label) so SD, VR, and Reconstruction can be told apart.
      const double bytesSent =
          static_cast<double>(cudaAwareCache.totalSendCount) * numEntries * sizeof(double);
      const double bytesRecv =
          static_cast<double>(cudaAwareCache.totalRecvCount) * numEntries * sizeof(double);

      pumipic::RecordTime(label + "_MPI_BytesSent_" + std::to_string(self), bytesSent);
      pumipic::RecordTime(label + "_MPI_BytesRecv_" + std::to_string(self), bytesRecv);

      // Post both the Irecv and the matching Isend for a proc together, in
      // the same loop iteration and under their own counts (recvCounts for
      // Irecv, sendCounts for Isend).
      //
      // The send is posted directly against meshField using that proc's
      // derived datatype (no sendDataGPU pointer/offset), unconditionally,
      // and recv requests are tracked separately from send requests (with a
      // parallel recvReqProc[] telling us which proc each recv request
      // belongs to) so the recvs can be drained incrementally with
      // MPI_Waitany below instead of all being blocked on together.
      int numNeighbors = 0;
      for(int proc = 0; proc < numProcsTot; proc++){
        if(proc == self) continue;
        bool hasComm = false;

        if(cudaAwareCache.recvCounts[proc] > 0){
          MPI_Request reqData;

          double* recvPtr =
              cudaAwareCache.recvDataGPU.data() +
              static_cast<size_t>(cudaAwareCache.recvOffsets[proc]) *
                  numEntries;

          mpiError = MPI_Irecv(
              recvPtr,
              cudaAwareCache.recvCounts[proc] * numEntries,
              MPI_DOUBLE,
              proc,
              2,
              comm,
              &reqData);
          if(mpiError != MPI_SUCCESS) break;
          recvRequests.push_back(reqData);
          recvReqProc.push_back(proc);
          hasComm = true;
        }

        if(cudaAwareCache.sendCounts[proc] > 0){
          MPI_Request reqData;

          mpiError = MPI_Isend(
              meshField.data(),
              1,
              cudaAwareCache.sendProcTypes[proc],
              proc,
              2,
              comm,
              &reqData);
          if(mpiError != MPI_SUCCESS) break;
          sendRequests.push_back(reqData);
          hasComm = true;
        }
         if(hasComm) numNeighbors++;
      }
      pumipic::RecordTime(label + "_MPI_NumNeighbors_" + std::to_string(self), static_cast<double>(numNeighbors));

      postTime = timer.seconds();

      pumipic::RecordTime(label + "_MPI_Post_" + std::to_string(self), postTime);

      timer.reset();

      // ---- Drain recvs with MPI_Waitany, unpacking each neighbor's
      // contribution the moment it lands instead of waiting for every
      // neighbor before unpacking any of them. Each unpack kernel launch is
      // asynchronous (no fence in the loop), so while neighbor i's data is
      // being scattered on the GPU, the CPU is already back inside
      // MPI_Waitany waiting on the rest - hiding unpack behind the wait for
      // stragglers rather than paying for it serially afterward. ----
      if(mpiError == MPI_SUCCESS && !recvRequests.empty()){
        auto recvIDGPUView = cudaAwareCache.recvIDGPU.view();
        auto recvDataGPUView = cudaAwareCache.recvDataGPU.view();

        for(size_t reqIdx = 0; reqIdx < recvRequests.size(); reqIdx++){
          int idx = MPI_UNDEFINED;

          mpiError = MPI_Waitany(
              static_cast<int>(recvRequests.size()),
              recvRequests.data(),
              &idx,
              MPI_STATUS_IGNORE);

          if(mpiError != MPI_SUCCESS || idx == MPI_UNDEFINED) break;

          const int proc = recvReqProc[idx];
          const int base = cudaAwareCache.recvOffsets[proc];
          const int count = cudaAwareCache.recvCounts[proc];

          if(op == 0){
            Kokkos::parallel_for(
                "halo add cached cuda-aware mpi incremental",
                Kokkos::RangePolicy<>(base, base + count),
                KOKKOS_LAMBDA(const int i){
                  const int vertex = recvIDGPUView(i);

                  for(int k = 0; k < numEntries; k++){
#ifdef POLYMPO_ASSUME_UNIQUE_HALO_CONTRIBS
                    meshField(vertex, k) +=
                        recvDataGPUView(i * numEntries + k);
#else
                    Kokkos::atomic_add(
                        &meshField(vertex, k),
                        recvDataGPUView(i * numEntries + k));
#endif
                  }
                });
          }
          else{
            Kokkos::parallel_for(
                "halo assign cached cuda-aware mpi incremental",
                Kokkos::RangePolicy<>(base, base + count),
                KOKKOS_LAMBDA(const int i){
                  const int vertex = recvIDGPUView(i);

                  for(int k = 0; k < numEntries; k++){
                    meshField(vertex, k) =
                        recvDataGPUView(i * numEntries + k);
                  }
                });
          }
          // Deliberately no fence here - see comment above the loop.
        }
      }

      // Sends don't feed any unpack step, but we still need to know they
      // completed before treating this call as done, so wait on them once,
      // after the recv/unpack loop rather than before it (so posting the
      // sends can't itself delay draining the recvs).
      if(mpiError == MPI_SUCCESS && !sendRequests.empty()){
        mpiError = MPI_Waitall(
            static_cast<int>(sendRequests.size()),
            sendRequests.data(),
            MPI_STATUSES_IGNORE);
      }

      waitTime = timer.seconds();

      pumipic::RecordTime(label + "_MPI_Waitall_" + std::to_string(self), waitTime);

      if(mpiError != MPI_SUCCESS){
        cudaAwareMPIDisabled = true;

        if(self == 0){
          std::cout
              << "[CUDA_AWARE_MPI] Batched device-pointer MPI failed."
              << std::endl;
        }

        if(recvRequests.empty() && sendRequests.empty()){
          if(self == 0){
            std::cout
                << "[CUDA_AWARE_MPI] Falling back to CPU-staged communication."
                << std::endl;
          }

          communicate_and_take_halo_contributions1(
              meshField,
              nEntities,
              numEntries,
              mode,
              op);
          return;
        }

        if(self == 0){
          std::cout
              << "[CUDA_AWARE_MPI] Failure happened after MPI requests were posted. "
              << "Set POLYMPO_DISABLE_CUDA_AWARE_MPI=1 before running to force the CPU-staged path."
              << std::endl;
        }

        MPI_Abort(comm, mpiError);
        return;
      }

      timer.reset();

      // Final fence: guarantees every unpack kernel launched inside the
      // MPI_Waitany loop above has actually finished before meshField is
      // used downstream. When neighbor arrivals are staggered, most of that
      // unpack work already finished while we were still waiting on
      // stragglers, so in that case this fence's cost is close to just the
      // last-arriving neighbor's unpack kernel rather than the sum of all
      // of them. If arrivals are tightly bunched instead, expect this to
      // look a lot like the old post-Waitall unpack cost.
      Kokkos::fence();

      if(mpiDiagnostics){
          pumipic::RecordTime(label + "_MPI_Diagnostics_Contribution_m" + std::to_string(mode) + "_e" +
          std::to_string(numEntries) + "_rank" + std::to_string(self), timer.seconds());
    }


    }

#else

    // Fallback path:
    // if CUDA_AWARE_MPI is not defined, use the original GPU-CPU staging function.
    template <typename ViewType>
    void communicate_and_take_halo_contributions1_improved(
        const ViewType& meshField,
        int nEntities,
        int numEntries,
        int mode,
        int op,
        const std::string& label){
      (void)label; // no per-call diagnostics on the CPU-staged fallback path

      communicate_and_take_halo_contributions1(
          meshField,
          nEntities,
          numEntries,
          mode,
          op);
    }

#endif

};

}//namespace polyMPO end

#endif
