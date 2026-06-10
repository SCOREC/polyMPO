#ifndef POLYMPO_MPM_H
#define POLYMPO_MPM_H

#include "pmpo_utils.hpp"
#include "pmpo_mesh.hpp"
#include "pmpo_materialPoints.hpp"

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
    std::vector<int>haloOwnerProcs;
    std::vector<std::vector<int>> haloOwnerLocalIDs;
    std::vector<std::vector<int>> ownerOwnerLocalIDs;
    std::vector<std::vector<int>> ownerHaloLocalIDs;

    void startCommunication();

    void communicate_and_take_halo_contributions(const Kokkos::View<double**>& meshField, int nEntities, int numEntries, int mode, int op);
    //Now Kokkos views are made 1D
    template <typename ViewType>
    void communicate_and_take_halo_contributions1(
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
      Kokkos::fence();
      pumipic::RecordTime("Communication-GPU to CPU-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());
      
      timer.reset();
      std::vector<double> fieldData1(nEntities*numEntries);
      for (int i=0;i<nEntities;i++)
        for (int j=0;j<numEntries;j++)
          fieldData1[i*numEntries+j]=reconVals_host(i,j);
      //std::memcpy(fieldData1.data(), reconVals_host.data(), nEntities*numEntries*sizeof(double));
      pumipic::RecordTime("New allocation + copy" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());




for (int i=0; i<nEntities; i++) {
  for (int j=0; j<numEntries; j++) {
    if (fieldData1[i*numEntries+j] != reconVals_host(i,j)) {
      std::cout
        << "Mismatch at "
        << i << " " << j
        << " : "
        << fieldData1[i*numEntries+j]
        << " "
        << reconVals_host(i,j)
        << std::endl;
      break;
    }
  }
}



      timer.reset();
      std::vector<std::vector<int>>    recvIDVec;
      std::vector<std::vector<double>> recvDataVec;
      communicateFields1(fieldData1, nEntities, numEntries, mode, recvIDVec, recvDataVec);
      pumipic::RecordTime("Communication-InterProcess-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());

      timer.reset();
      int numProcsTot =  recvIDVec.size();
      //Flatten IDs 
      int totalSize = 0;
      std::vector<int> offsets(numProcsTot, 0); 
      for(int i=0; i<numProcsTot; i++) {
        offsets[i] = totalSize;
        totalSize += recvIDVec[i].size();
      }
      std::vector<int> flatIDVec(totalSize, 0);
      for(int i=0; i<numProcsTot; i++) {
        std::copy(recvIDVec[i].begin(), recvIDVec[i].end(), flatIDVec.begin() + offsets[i]);
      }

      Kokkos::View<int*> recvIDGPU("recvIDGPU", totalSize);
      auto hostView = Kokkos::View<int*, Kokkos::HostSpace>("recvIDCPU", totalSize);
      std::copy(flatIDVec.begin(), flatIDVec.end(), hostView.data());
      Kokkos::deep_copy(recvIDGPU, hostView);

      //Flatten Data
      int totalSize_data=0;
      std::vector<int> offsets_data(numProcsTot, 0);
      for(int i=0; i<numProcsTot; i++) {
        offsets_data[i] = totalSize_data;
        totalSize_data += recvDataVec[i].size();
      }
      std::vector<double> flatDataVec(totalSize_data, 0);
      for(int i=0; i<numProcsTot; i++) {
        std::copy(recvDataVec[i].begin(), recvDataVec[i].end(), flatDataVec.begin() + offsets_data[i]);
      }
      Kokkos::View<double*> recvDataGPU("recvDataGPU", totalSize_data);
      auto hostView_data= Kokkos::View<double*, Kokkos::HostSpace>("recvDataCPU", totalSize_data);
      std::copy(flatDataVec.begin(), flatDataVec.end(), hostView_data.data()); 
      Kokkos::deep_copy(recvDataGPU, hostView_data);
      Kokkos::fence();
      //Assertions
      assert(totalSize_data == totalSize*numEntries);
      for (int i=0; i<numProcsTot; i++){
        assert(recvDataVec[i].size() == recvIDVec[i].size() * numEntries);
      }
      pumipic::RecordTime("Communication-CPU to GPU-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());
  
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
      pumipic::RecordTime("Communication-GPU reduction-E-" + std::to_string(numEntries) + "-" + std::to_string(self), timer.seconds());
    }


    void communicateFields(const std::vector<std::vector<double>>& fieldData, const int numEntities, const int numEntries, int mode,
                              std::vector<std::vector<int>>& recvIDVec,  std::vector<std::vector<double>>& recvDataVec);
    //1D Kokkos view will be copied direcctly to 1D std::vector
    void communicateFields1(
        const std::vector<double>& fieldData, 
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
            sendDataVec[ownerProc].push_back(fieldData[(numOwnersTot+iEnt)*numEntries + iDouble]);
        }
      }
      else if(mode == 1){
        // Owner sends to halos
        for (size_t iProc=0; iProc<ownerOwnerLocalIDs.size(); iProc++) {
          for (auto& ownerID : ownerOwnerLocalIDs[iProc]) {
            for (int iDouble = 0; iDouble < numEntries; iDouble++)
              sendDataVec[iProc].push_back(fieldData[ownerID*numEntries + iDouble]);
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
};

}//namespace polyMPO end

#endif

