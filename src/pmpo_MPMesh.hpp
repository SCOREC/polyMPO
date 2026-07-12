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
    // via .view(). Used only for the 4 buffers below that get handed
    // directly to MPI_Isend/Irecv under CUDA-aware MPI.
    //
    // Why not just a Kokkos::View<T*, Kokkos::CudaSpace>: when Kokkos is
    // built with Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=ON (the default since
    // Kokkos 4.2), View allocations use cudaMallocAsync/memory pools.
    // cuIpcGetMemHandle (which Cray MPICH/GTL uses for intra-node GPU-to-GPU
    // sends) rejects pool allocations with CUDA_ERROR_INVALID_VALUE. Since
    // polyMPO isn't allowed to touch the Kokkos build config, these 4
    // buffers bypass Kokkos's allocator entirely via a direct cudaMalloc,
    // which cuIpcGetMemHandle always accepts, regardless of how the rest of
    // Kokkos (or the rest of the app's Views) is configured.
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
    // Every neighbor's data lives in one shared allocation (see
    // CudaAwareMPIFieldCache below) and MPI is given "base pointer + byte
    // offset" per neighbor rather than a separate allocation per neighbor.
    // If you ever need to fall back to one allocation per neighbor (e.g. an
    // MPI/GPU stack that mishandles offset device pointers for CUDA IPC),
    // restore the per-proc-buffer version from version control.
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

      // Single batched GPU buffers (one allocation each, instead of one
      // Kokkos::View per neighbor proc). Per-proc slices are
      // [offset, offset + count) for the ID buffers, and
      // [offset * numEntries, (offset + count) * numEntries) for the data
      // buffers. MPI is given "buffer base pointer + offset", not a
      // separate allocation per proc.
      CudaAwareMPIIntBuffer sendEntityGPU;
      CudaAwareMPIIntBuffer recvIDGPU;

      CudaAwareMPIDoubleBuffer sendDataGPU;
      CudaAwareMPIDoubleBuffer recvDataGPU;
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

    // Fully CUDA-aware MPI version, batched buffer variant:
    // Field data is sent/received using GPU pointers. Receive IDs are cached
    // once from the fixed halo/owner mapping and are not sent every call.
    //
    // Every neighbor's send/recv entity-ID list and data live in ONE big
    // GPU buffer each (laid out back-to-back in proc order), instead of one
    // Kokkos::View allocation per neighbor. Packing/unpacking is a single
    // kernel launch over all neighbors' entities at once instead of one
    // launch per neighbor, and MPI_Isend/Irecv use "buffer base pointer +
    // offset" into that single buffer per proc. This is what actually
    // shrinks MPI_Wait time: fewer, larger, more uniform in-flight
    // transfers instead of many small independent ones.
    //
    // Note: an earlier version of this cache used one Kokkos::View
    // allocation per neighbor specifically to avoid handing MPI a
    // "base pointer + offset" GPU address, out of concern for CUDA IPC
    // issues on Cray MPICH/GTL. That failure mode was root-caused to
    // Kokkos allocating device Views via cudaMallocAsync (invalid for
    // cuIpcGetMemHandle), not to offset pointers themselves, and is fixed
    // by building Kokkos with -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF.
    // If you ever do hit IPC trouble that tracks back to offset pointers
    // specifically, the per-proc-buffer version can be restored from
    // version control.
    //
    // Important:
    // This function caches communication metadata and GPU buffers per
    // (mode, numEntries). If the communication pattern changes, clear
    // cudaAwareMPICaches before the next call.
    template <typename ViewType>
    void communicate_and_take_halo_contributions1_improved(
        const ViewType& meshField,
        int nEntities,
        int numEntries,
        int mode,
        int op){

      int self, numProcsTot;

      MPI_Comm comm = p_MPs->getMPIComm();

      MPI_Comm_rank(comm, &self);
      MPI_Comm_size(comm, &numProcsTot);

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
            << "[CUDA_AWARE_MPI] Using batched single-buffer GPU-aware MPI path in communicate_and_take_halo_contributions1_improved()"
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

        cudaAwareCache.sendEntityGPU.allocate(totalSend);
        cudaAwareCache.sendDataGPU.allocate(totalSend * numEntries);
        cudaAwareCache.recvIDGPU.allocate(totalRecv);
        cudaAwareCache.recvDataGPU.allocate(totalRecv * numEntries);

        // ---- Build the flattened send-entity list (host, then one deep_copy) ----
        if(totalSend > 0){
          auto sendEntityCPU =
              Kokkos::View<int*, Kokkos::HostSpace>(
                  "sendEntityCPU_batched", totalSend);

          if(mode == 0){
            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;
              if(cudaAwareCache.sendCounts[proc] <= 0) continue;

              assert(haloOwnerLocalIDs[proc].size() ==
                     static_cast<size_t>(cudaAwareCache.sendCounts[proc]));
            }

            std::vector<int> cursor(cudaAwareCache.sendOffsets);

            for(int iEnt = 0; iEnt < numHalosTot; iEnt++){
              int ownerProc = haloOwnerProcs[iEnt];
              if(ownerProc == self) continue;

              sendEntityCPU(cursor[ownerProc]) = numOwnersTot + iEnt;
              cursor[ownerProc]++;
            }

            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;

              assert(cursor[proc] ==
                     cudaAwareCache.sendOffsets[proc] +
                         cudaAwareCache.sendCounts[proc]);
            }
          }
          else{
            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;

              int sendCount = cudaAwareCache.sendCounts[proc];
              if(sendCount <= 0) continue;

              assert(ownerOwnerLocalIDs[proc].size() ==
                     static_cast<size_t>(sendCount));

              int base = cudaAwareCache.sendOffsets[proc];

              for(int i = 0; i < sendCount; i++){
                sendEntityCPU(base + i) = ownerOwnerLocalIDs[proc][i];
              }
            }
          }

          Kokkos::deep_copy(cudaAwareCache.sendEntityGPU.view(), sendEntityCPU);
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

        pumipic::RecordTime(
            "SD: CUDA-aware MPI Cache Build m" + std::to_string(mode) +
                " e" + std::to_string(numEntries) + "-" + std::to_string(self),
            timer.seconds());

        timer.reset();
      }

      // ---- Pack: ONE kernel over all neighbors' send entities at once ----
      if(cudaAwareCache.totalSendCount > 0){
        auto sendEntityGPU = cudaAwareCache.sendEntityGPU.view();
        auto sendDataGPU = cudaAwareCache.sendDataGPU.view();

        Kokkos::parallel_for(
            "pack cached cuda-aware mpi send buffer batched",
            cudaAwareCache.totalSendCount,
            KOKKOS_LAMBDA(const int i){
              int entity = sendEntityGPU(i);

              for(int k = 0; k < numEntries; k++){
                sendDataGPU(i * numEntries + k) =
                    meshField(entity, k);
              }
            });
      }

      Kokkos::fence();

      pumipic::RecordTime(
          "SD: CUDA-aware MPI Pack m" + std::to_string(mode) +
              " e" + std::to_string(numEntries) + "-" + std::to_string(self),
          timer.seconds());

      timer.reset();

      Kokkos::Timer mpiTotalTimer;
      std::vector<MPI_Request> requests;
      requests.reserve(2 * numProcsTot);
      int mpiError = MPI_SUCCESS;

      for(int proc = 0; proc < numProcsTot; proc++){
        if(proc == self) continue;

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
          requests.push_back(reqData);
        }

        if(cudaAwareCache.sendCounts[proc] > 0){
          MPI_Request reqData;

          double* sendPtr =
              cudaAwareCache.sendDataGPU.data() +
              static_cast<size_t>(cudaAwareCache.sendOffsets[proc]) *
                  numEntries;

          mpiError = MPI_Isend(
              sendPtr,
              cudaAwareCache.sendCounts[proc] * numEntries,
              MPI_DOUBLE,
              proc,
              2,
              comm,
              &reqData);
          if(mpiError != MPI_SUCCESS) break;
          requests.push_back(reqData);
        }
      }

      pumipic::RecordTime(
          "SD: CUDA-aware MPI Post m" + std::to_string(mode) +
              " e" + std::to_string(numEntries) + "-" + std::to_string(self),
          timer.seconds());

      timer.reset();

      if(mpiError == MPI_SUCCESS && !requests.empty()){
        mpiError = MPI_Waitall(
            static_cast<int>(requests.size()),
            requests.data(),
            MPI_STATUSES_IGNORE);
      }

      pumipic::RecordTime(
          "SD: CUDA-aware MPI Wait m" + std::to_string(mode) +
              " e" + std::to_string(numEntries) + "-" + std::to_string(self),
          timer.seconds());

      if(mpiError != MPI_SUCCESS){
        cudaAwareMPIDisabled = true;

        if(self == 0){
          std::cout
              << "[CUDA_AWARE_MPI] Batched device-pointer MPI failed."
              << std::endl;
        }

        if(requests.empty()){
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

        timer.reset();

        if(self == 0){
          std::cout
              << "[CUDA_AWARE_MPI] Failure happened after MPI requests were posted. "
              << "Set POLYMPO_DISABLE_CUDA_AWARE_MPI=1 before running to force the CPU-staged path."
              << std::endl;
        }

        MPI_Abort(comm, mpiError);
        return;
      }

      pumipic::RecordTime(
          "SD: CUDA-aware MPI Comm m" + std::to_string(mode) +
              " e" + std::to_string(numEntries) + "-" + std::to_string(self),
          mpiTotalTimer.seconds());

      timer.reset();

      // ---- Unpack: ONE kernel over all neighbors' recv entities at once ----
      if(cudaAwareCache.totalRecvCount > 0){
        auto recvIDGPU = cudaAwareCache.recvIDGPU.view();
        auto recvDataGPU = cudaAwareCache.recvDataGPU.view();

        if(op == 0){
          Kokkos::parallel_for(
              "halo add cached cuda-aware mpi batched",
              cudaAwareCache.totalRecvCount,
              KOKKOS_LAMBDA(const int i){
                const int vertex = recvIDGPU(i);

                for(int k = 0; k < numEntries; k++){
#ifdef POLYMPO_ASSUME_UNIQUE_HALO_CONTRIBS
                  meshField(vertex, k) +=
                      recvDataGPU(i * numEntries + k);
#else
                  Kokkos::atomic_add(
                      &meshField(vertex, k),
                      recvDataGPU(i * numEntries + k));
#endif
                }
              });
        }
        else{
          Kokkos::parallel_for(
              "halo assign cached cuda-aware mpi batched",
              cudaAwareCache.totalRecvCount,
              KOKKOS_LAMBDA(const int i){
                const int vertex = recvIDGPU(i);

                for(int k = 0; k < numEntries; k++){
                  meshField(vertex, k) =
                      recvDataGPU(i * numEntries + k);
                }
              });
        }
      }

      Kokkos::fence();

      pumipic::RecordTime(
          "SD: CUDA-aware MPI Contribution m" + std::to_string(mode) +
              " e" + std::to_string(numEntries) + "-" + std::to_string(self),
          timer.seconds());
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
        int op){

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
