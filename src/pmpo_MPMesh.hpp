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
    void communicate_and_take_halo_contributions(const Kokkos::View<double**>& meshField, int nEntities, int numEntries, int mode, int op);

    template <typename ViewType>
    void communicate_and_take_halo_contributions_staged(
        const ViewType& meshField,
        int nEntities,
        int numEntries,
        int mode ,
        int op){

      int self;
      MPI_Comm comm = p_MPs->getMPIComm();
      MPI_Comm_rank(comm, &self);

      Kokkos::Timer timer;
      auto reconVals_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), meshField);
      pumipic::RecordTime("SD: GPU-CPU copy-" + std::to_string(self), timer.seconds());

      timer.reset();
      std::vector<std::vector<int>>    recvIDVec;
      std::vector<std::vector<double>> recvDataVec;
      pumipic::RecordTime("SD: Recv Vec Allocation-" + std::to_string(self), timer.seconds());
      timer.reset();
      //communicateFieldsFromHostView(fieldData1, nEntities, numEntries, mode, recvIDVec, recvDataVec);
      communicateFieldsFromHostView(reconVals_host, nEntities, numEntries, mode, recvIDVec, recvDataVec);
      pumipic::RecordTime("SD: IP Comm-" + std::to_string(self), timer.seconds());

      timer.reset();
      int numProcsTot = recvIDVec.size();
      //Flatten IDs
      int totalSize = 0;
      std::vector<int> offsets(numProcsTot, 0); 
      for(int i=0; i<numProcsTot; i++) {
        offsets[i] = totalSize;
        totalSize += recvIDVec[i].size();
      }
      std::vector<int> flatIDVec(totalSize, 0);
      for(int i=0; i<numProcsTot; i++){
        std::copy(recvIDVec[i].begin(), recvIDVec[i].end(), flatIDVec.begin() + offsets[i]);
      }
      pumipic::RecordTime("SD: Flatten IDs-" + std::to_string(self), timer.seconds());

      timer.reset();
      Kokkos::View<int*> recvIDGPU("recvIDGPU", totalSize);
      auto hostView = Kokkos::View<int*, Kokkos::HostSpace>("recvIDCPU", totalSize);
      std::copy(flatIDVec.begin(), flatIDVec.end(), hostView.data());
      Kokkos::deep_copy(recvIDGPU, hostView);
      Kokkos::fence();
      pumipic::RecordTime("SD: Copy CPU-GPU-" + std::to_string(self), timer.seconds());

      //Flatten Data
      timer.reset();
      int totalSize_data=0;
      std::vector<int> offsets_data(numProcsTot, 0);
      for(int i=0; i<numProcsTot; i++){
        offsets_data[i] = totalSize_data;
        totalSize_data += recvDataVec[i].size();
      }
      std::vector<double> flatDataVec(totalSize_data, 0);
      for(int i=0; i<numProcsTot; i++) {
        std::copy(recvDataVec[i].begin(), recvDataVec[i].end(), flatDataVec.begin() + offsets_data[i]);
      }
      pumipic::RecordTime("SD: Flatten Data-" + std::to_string(self), timer.seconds());

      timer.reset();
      Kokkos::View<double*> recvDataGPU("recvDataGPU", totalSize_data);
      auto hostView_data= Kokkos::View<double*, Kokkos::HostSpace>("recvDataCPU", totalSize_data);
      std::copy(flatDataVec.begin(), flatDataVec.end(), hostView_data.data()); 
      Kokkos::deep_copy(recvDataGPU, hostView_data);
      Kokkos::fence();
      assert(totalSize_data == totalSize*numEntries);
      for (int i=0; i<numProcsTot; i++){
        assert(recvDataVec[i].size() == recvIDVec[i].size() * numEntries);
      }
      pumipic::RecordTime("SD: Copy CPU-GPU2-" + std::to_string(self), timer.seconds());
      //Take contributions from other procs
      timer.reset();
      Kokkos::parallel_for("halo contribution", recvIDGPU.size(), KOKKOS_LAMBDA(const int i){
        int vertex = recvIDGPU(i);
        for(int k=0; k<numEntries; k++){
          if(op==0) Kokkos::atomic_add(&meshField(vertex,k), recvDataGPU(i*numEntries+k));
          if(op==1) meshField(vertex, k) = recvDataGPU(i * numEntries + k);
        }
      });
      Kokkos::fence();
      pumipic::RecordTime("SD: Contribution" + std::to_string(self), timer.seconds());
    }

    void communicateFields(const std::vector<std::vector<double>>& fieldData, const int numEntities, const int numEntries, int mode,
                              std::vector<std::vector<int>>& recvIDVec, std::vector<std::vector<double>>& recvDataVec);

    template <class ViewType>
    void communicateFieldsFromHostView(
        const ViewType& fieldData,
        const int numEntities, const int numEntries, int mode,
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
        if(i==self) continue;

        int numToSend = 0, numToRecv = 0; 
        if(mode == 0) {
          //gather (halos send to owners)
          numToSend = numOwnersOnOtherProcs[i];
          numToRecv = numHalosOnOtherProcs[i];
        }
        else{ 
          //scatter (owners send to halos)
          numToSend = numHalosOnOtherProcs[i];
          numToRecv = numOwnersOnOtherProcs[i];
        }

        if(numToSend > 0){
          sendDataVec[i].reserve(numToSend*numEntries);
        }
        if(numToRecv > 0){
          recvDataVec[i].resize(numToRecv*numEntries);
          recvIDVec[i].resize(numToRecv);
        }
      }

      if(mode == 0){
        // Halos sends to owners
        for (int iEnt = 0; iEnt < numHalosTot; iEnt++){
          auto ownerProc = haloOwnerProcs[iEnt];
          for (int iDouble = 0; iDouble < numEntries; iDouble++)
            sendDataVec[ownerProc].push_back(fieldData(numOwnersTot+iEnt, iDouble));
        }
      }
      else if(mode == 1){
        // Owner sends to halos
        for (size_t iProc=0; iProc<ownerOwnerLocalIDs.size(); iProc++) {
          for (auto& ownerID : ownerOwnerLocalIDs[iProc]) {
            for (int iDouble = 0; iDouble < numEntries; iDouble++)
              sendDataVec[iProc].push_back(fieldData(ownerID, iDouble));
          }
        }
      }

      std::vector<MPI_Request> requests;
      requests.reserve(4*numProcsTot);
      for(int proc = 0; proc < numProcsTot; proc++){
        if(proc == self) continue; 
        if(mode == 0 && numHalosOnOtherProcs[proc]){
          assert(recvIDVec[proc].size() == (size_t)numHalosOnOtherProcs[proc]);
          assert(recvDataVec[proc].size() == recvIDVec[proc].size() * (size_t)numEntries);
          MPI_Request req3, req4;
          MPI_Irecv(recvIDVec[proc].data(), recvIDVec[proc].size(), MPI_INT, proc, 1, comm, &req3);
          MPI_Irecv(recvDataVec[proc].data(), recvDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req4);
          requests.push_back(req3);
          requests.push_back(req4);
        }
        if(mode == 0 && numOwnersOnOtherProcs[proc]) {
          assert(haloOwnerLocalIDs[proc].size() == (size_t)numOwnersOnOtherProcs[proc]);
          assert(sendDataVec[proc].size() == haloOwnerLocalIDs[proc].size() * (size_t)numEntries);
          MPI_Request req1, req2;
          MPI_Isend(haloOwnerLocalIDs[proc].data(), haloOwnerLocalIDs[proc].size(), MPI_INT, proc, 1, comm, &req1);
          MPI_Isend(sendDataVec[proc].data(), sendDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req2);
          requests.push_back(req1);
          requests.push_back(req2);
        }

        if(mode == 1 && numOwnersOnOtherProcs[proc]){
          MPI_Request req3, req4;
          MPI_Irecv(recvIDVec[proc].data(), recvIDVec[proc].size(), MPI_INT, proc, 1, comm, &req3);
          MPI_Irecv(recvDataVec[proc].data(), recvDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req4);
          requests.push_back(req3);
          requests.push_back(req4);
        }
        if(mode == 1 && numHalosOnOtherProcs[proc]) {
          MPI_Request req1, req2;
          MPI_Isend(ownerHaloLocalIDs[proc].data(), ownerHaloLocalIDs[proc].size(), MPI_INT, proc, 1, comm, &req1);
          MPI_Isend(sendDataVec[proc].data(), sendDataVec[proc].size(), MPI_DOUBLE, proc, 2, comm, &req2);
          requests.push_back(req1);
          requests.push_back(req2);
        }
      }
      MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    }
   
    MPMesh(Mesh* inMesh, MaterialPoints* inMPs):
      p_mesh(inMesh), p_MPs(inMPs) {
    };

    ~MPMesh() {
      delete p_mesh;
      delete p_MPs;
    }

    //MP advection and tracking
    void CVTTrackingEdgeCenterBased(Vec2dView dx);
    void CVTTrackingElmCenterBased(const int printVTPIndex = -1);
    void T2LTracking(Vec2dView dx);
    bool push1P();
    void push_ahead();
    void push_swap();
    void push_swap_pos();
    void push();

    //Used before advection to interpolate fields from mesh to MPs
    //And also before reconstruction
    void calcBasis();

    //Reconstruction
    DoubleView assemblyV0();
    template <MeshFieldIndex meshFieldIndex>
    void assemblyVtx0();
    template <MeshFieldIndex meshFieldIndex>
    void assemblyElm0();
    template <MeshFieldIndex meshFieldIndex>
    void assemblyVtx1();
    void reconstruct_coeff_full();
    void invertMatrix(const Kokkos::View<double**>& vtxMatrices, const double& radius);
    Kokkos::View<double*[vec3d_nEntries][vec4d_nEntries]> precomputedVtxCoeffs_new;
    Kokkos::View<double*> nearAnEdge;  
    Kokkos::View<double*> vtxMatrixMass;

    //Not used currently
    std::map<MeshFieldIndex, std::function<void()>> reconstructSlice = std::map<MeshFieldIndex, std::function<void()>>();
    template <MaterialPointSlice index>
    DoubleView wtScaAssembly();
    template <MaterialPointSlice index>
    Vec2dView wtVec2Assembly();
    template <MeshFieldIndex meshFieldIndex>
    void assembly(int order, MeshFieldType type, bool basisWeightFlag, bool massWeightFlag);
    template<MeshFieldIndex meshFieldIndex>
    void setReconstructSlice(int order, MeshFieldType type);
    void reconstructSlices();

    void printVTP_mesh(int printVTPIndex);
    void writeMPTrackingVTP(int printVTPIndex, int numMPs, const Vec3dView& history, const Vec3dView& resultLeft,
                            const Vec3dView& resultRight, const Vec3dView& mpTgtPosArray);
    void calculateStrain();
    void calculateStress(const int constitutive_relation);
    void calculateStressDivergence();

    //Prints the selected halo-exchange path once, on rank 0
    bool haloExchangePathReported = false;
    void reportHaloExchangePath(const char* path){
      if(haloExchangePathReported) return;
      haloExchangePathReported = true;
      int self;
      MPI_Comm_rank(p_MPs->getMPIComm(), &self);
      if(self == 0){
        std::cout << "polyMPO halo exchange: " << path << std::endl;
      }
    }

#ifdef GPU_AWARE_MPI

#ifdef KOKKOS_ENABLE_CUDA
    //Device buffer allocated with plain cudaMalloc and exposed as an
    //unmanaged Kokkos::View. Used for buffers passed directly to MPI.
    //
    //Kokkos (>= 4.2, CUDA_MALLOC_ASYNC=ON by default) allocates Views from
    //cudaMallocAsync memory pools, which cuIpcGetMemHandle rejects
    //(CUDA_ERROR_INVALID_VALUE). Cray MPICH/GTL relies on CUDA IPC for
    //intra-node GPU-to-GPU transfers, so these buffers bypass the Kokkos
    //allocator regardless of how Kokkos was built.
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
            throw std::runtime_error(std::string("RawCudaMPIBuffer: cudaMalloc failed: ") + cudaGetErrorString(err));
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
        return Kokkos::View<T*, Kokkos::CudaSpace, Kokkos::MemoryUnmanaged>(ptr, count);
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
    //Non-CUDA backend: plain Kokkos::View with the same interface as RawCudaMPIBuffer
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

    bool cudaAwareMPICacheValid = true;
    bool cudaAwareMPIDisabled = false;
    bool mpiGpuSupportChecked = false;
    bool mpiGpuSupportEnabled = false;

    struct CudaAwareMPIFieldCache{
      bool valid = false;
      int cachedNumProcs = -1;

      std::vector<int> sendCounts;
      std::vector<int> recvCounts;
      std::vector<int> sendOffsets; //prefix sum of sendCounts, in entities
      std::vector<int> recvOffsets; //prefix sum of recvCounts, in entities

      int totalSendCount = 0;
      int totalRecvCount = 0;

      //Batched buffers for all neighbors. Per-proc slices are
      //[offset, offset + count) for IDs and
      //[offset * numEntries, (offset + count) * numEntries) for data.
      CudaAwareMPIIntBuffer sendEntityGPU;
      CudaAwareMPIIntBuffer recvIDGPU;
      CudaAwareMPIDoubleBuffer sendDataGPU;
      CudaAwareMPIDoubleBuffer recvDataGPU;
    };

    std::map<std::pair<int, int>, CudaAwareMPIFieldCache> cudaAwareMPICaches;

    //The halo-exchange path follows MPICH_GPU_SUPPORT_ENABLED:
    //  1                                  -> GPU-aware MPI
    //  anything else (0, unset, yes, ...) -> CPU-staged path
    //MPI must not be given device pointers unless GPU support is enabled;
    //doing so crashes instead of returning an MPI error.
    bool mpiGpuSupport(){
      if(!mpiGpuSupportChecked){
        const char* value = std::getenv("MPICH_GPU_SUPPORT_ENABLED");
        mpiGpuSupportEnabled = value != nullptr && std::string(value) == "1";
        mpiGpuSupportChecked = true;
      }
      return mpiGpuSupportEnabled;
    }

    //GPU-aware MPI path: field data is sent/received directly from GPU
    //buffers. Receive IDs come from the fixed halo/owner mapping and are
    //not exchanged. All neighbors share one batched buffer per direction,
    //packed/unpacked with a single kernel each.
    //
    //Metadata and buffers are cached per (mode, numEntries). If the
    //communication pattern changes, clear cudaAwareMPICaches first.
    template <typename ViewType>
    void communicate_and_take_halo_contributions_gpu_aware(
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

      if(cudaAwareMPIDisabled || !mpiGpuSupport()){
        reportHaloExchangePath("CPU-staged (MPICH_GPU_SUPPORT_ENABLED is not 1)");
        communicate_and_take_halo_contributions_staged(meshField, nEntities, numEntries, mode, op);
        return;
      }

      reportHaloExchangePath("GPU-aware MPI");

      if(!cudaAwareMPICacheValid){
        cudaAwareMPICaches.clear();
        cudaAwareMPICacheValid = true;
      }

      auto& cudaAwareCache = cudaAwareMPICaches[std::make_pair(mode, numEntries)];
      const bool needRebuild = (!cudaAwareCache.valid) || (cudaAwareCache.cachedNumProcs != numProcsTot);

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

        //Build flattened send-entity list on host, then copy to device
        if(totalSend > 0){
          auto sendEntityCPU = Kokkos::View<int*, Kokkos::HostSpace>("sendEntityCPU_batched", totalSend);
          if(mode == 0){
            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;
              if(cudaAwareCache.sendCounts[proc] <= 0) continue;
              assert(haloOwnerLocalIDs[proc].size() == (size_t)cudaAwareCache.sendCounts[proc]);
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
              assert(cursor[proc] == cudaAwareCache.sendOffsets[proc] + cudaAwareCache.sendCounts[proc]);
            }
          }
          else{
            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;
              int sendCount = cudaAwareCache.sendCounts[proc];
              if(sendCount <= 0) continue;
              assert(ownerOwnerLocalIDs[proc].size() == (size_t)sendCount);
              int base = cudaAwareCache.sendOffsets[proc];
              for(int i = 0; i < sendCount; i++){
                sendEntityCPU(base + i) = ownerOwnerLocalIDs[proc][i];
              }
            }
          }
          Kokkos::deep_copy(cudaAwareCache.sendEntityGPU.view(), sendEntityCPU);
        }

        //Build flattened recv-ID list on host, then copy to device
        if(totalRecv > 0){
          auto recvIDCPU = Kokkos::View<int*, Kokkos::HostSpace>("recvIDCPU_batched", totalRecv);
          if(mode == 0){
            for(int proc = 0; proc < numProcsTot; proc++){
              if(proc == self) continue;
              int recvCount = cudaAwareCache.recvCounts[proc];
              if(recvCount <= 0) continue;
              assert(ownerOwnerLocalIDs[proc].size() == (size_t)recvCount);
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
              assert(cursor[proc] == cudaAwareCache.recvOffsets[proc] + cudaAwareCache.recvCounts[proc]);
            }
          }
          Kokkos::deep_copy(cudaAwareCache.recvIDGPU.view(), recvIDCPU);
        }

        cudaAwareCache.valid = true;
      }

      //Pack send buffer
      if(cudaAwareCache.totalSendCount > 0){
        auto sendEntityGPU = cudaAwareCache.sendEntityGPU.view();
        auto sendDataGPU = cudaAwareCache.sendDataGPU.view();
        Kokkos::parallel_for("pack cached gpu-aware mpi send buffer batched", cudaAwareCache.totalSendCount, KOKKOS_LAMBDA(const int i){
          int entity = sendEntityGPU(i);
          for(int k = 0; k < numEntries; k++){
            sendDataGPU(i * numEntries + k) = meshField(entity, k);
          }
        });
      }
      //MPI is not stream-aware: packing must finish before Isend
      Kokkos::fence();

      std::vector<MPI_Request> requests;
      requests.reserve(2 * numProcsTot);
      int mpiError = MPI_SUCCESS;

      //Post Irecv/Isend per neighbor
      for(int proc = 0; proc < numProcsTot; proc++){
        if(proc == self) continue;

        if(cudaAwareCache.recvCounts[proc] > 0){
          MPI_Request reqData;
          double* recvPtr = cudaAwareCache.recvDataGPU.data() + (size_t)cudaAwareCache.recvOffsets[proc] * numEntries;
          mpiError = MPI_Irecv(recvPtr, cudaAwareCache.recvCounts[proc] * numEntries, MPI_DOUBLE, proc, 2, comm, &reqData);
          if(mpiError != MPI_SUCCESS) break;
          requests.push_back(reqData);
        }

        if(cudaAwareCache.sendCounts[proc] > 0){
          MPI_Request reqData;
          double* sendPtr = cudaAwareCache.sendDataGPU.data() + (size_t)cudaAwareCache.sendOffsets[proc] * numEntries;
          mpiError = MPI_Isend(sendPtr, cudaAwareCache.sendCounts[proc] * numEntries, MPI_DOUBLE, proc, 2, comm, &reqData);
          if(mpiError != MPI_SUCCESS) break;
          requests.push_back(reqData);
        }
      }

      if(mpiError == MPI_SUCCESS && !requests.empty()){
        mpiError = MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      }

      if(mpiError != MPI_SUCCESS){
        cudaAwareMPIDisabled = true;
        if(self == 0){
          std::cout << "[GPU_AWARE_MPI] Batched device-pointer MPI failed." << std::endl;
        }

        if(requests.empty()){
          if(self == 0){
            std::cout << "[GPU_AWARE_MPI] Falling back to CPU-staged communication." << std::endl;
          }
          communicate_and_take_halo_contributions_staged(meshField, nEntities, numEntries, mode, op);
          return;
        }

        if(self == 0){
          std::cout << "[GPU_AWARE_MPI] Failure happened after MPI requests were posted. "
                    << "Set MPICH_GPU_SUPPORT_ENABLED=0 before running to use the CPU-staged path." << std::endl;
        }
        MPI_Abort(comm, mpiError);
        return;
      }

      //Unpack received contributions
      if(cudaAwareCache.totalRecvCount > 0){
        auto recvIDGPU = cudaAwareCache.recvIDGPU.view();
        auto recvDataGPU = cudaAwareCache.recvDataGPU.view();
        if(op == 0){
          Kokkos::parallel_for("halo add cached gpu-aware mpi batched", cudaAwareCache.totalRecvCount, KOKKOS_LAMBDA(const int i){
            const int vertex = recvIDGPU(i);
            for(int k = 0; k < numEntries; k++){
#ifdef POLYMPO_ASSUME_UNIQUE_HALO_CONTRIBS
              meshField(vertex, k) += recvDataGPU(i * numEntries + k);
#else
              Kokkos::atomic_add(&meshField(vertex, k), recvDataGPU(i * numEntries + k));
#endif
            }
          });
        }
        else{
          Kokkos::parallel_for("halo assign cached gpu-aware mpi batched", cudaAwareCache.totalRecvCount, KOKKOS_LAMBDA(const int i){
            const int vertex = recvIDGPU(i);
            for(int k = 0; k < numEntries; k++){
              meshField(vertex, k) = recvDataGPU(i * numEntries + k);
            }
          });
        }
      }
      Kokkos::fence();
    }

#else

    //GPU_AWARE_MPI not defined: use the CPU-staged path
    template <typename ViewType>
    void communicate_and_take_halo_contributions_gpu_aware(
        const ViewType& meshField,
        int nEntities,
        int numEntries,
        int mode,
        int op){

      reportHaloExchangePath("CPU-staged (built without GPU_AWARE_MPI)");
      communicate_and_take_halo_contributions_staged(meshField, nEntities, numEntries, mode, op);
    }

#endif

};

}//namespace polyMPO end

#endif
