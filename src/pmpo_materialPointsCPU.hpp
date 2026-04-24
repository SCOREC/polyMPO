#include "pmpo_mesh.hpp"
#include "pmpo_materialPoints.hpp"

namespace polyMPO{

enum CPU_MP_Slice{
  MPF_IceAreaMP = 0,
  MPF_IceEnthalpy,
  MPF_NUM_CPU_FIELDS
};

class MaterialPointsCPU{
  
  private:
    MaterialPoints* mps_g;

  public:

    MaterialPointsCPU()=default;

    MaterialPointsCPU(MaterialPoints* p_MPs):mps_g(p_MPs){};

    ~MaterialPointsCPU()=default;


    std::vector<std::vector<double>> CPU_fields;

    void setFields(const int nParticles, const double* field_array){  
      
      int nRows = CPU_fields.size();
      if (nRows != nParticles){
         CPU_fields.resize(nParticles, std::vector<double>(MPF_NUM_CPU_FIELDS));
      }
      for (auto mp = 0; mp < nRows; ++mp)
        CPU_fields[mp][0] = field_array[mp];
    }

    void migrateParticlesCPU(){
      
      int numMPs = mps_g->getCount();
      auto mpAppID  = mps_g->getData<MPF_MP_APP_ID>();
      auto MPs2Proc = mps_g->getData<MPF_Tgt_Proc_ID>();
      
      Kokkos::View<int*> mpAppID_array("mpAppID", mps_g->getCapacity());
      Kokkos::View<int*> mpProc_array("mpProc", mps_g->getCapacity());
      Kokkos::View<int*> counter("counter", 1);
      Kokkos::deep_copy(counter, 0); 
      // Store MPAppID and Process as continuous Arrays
      auto setMPProcAppID = PS_LAMBDA(const int& e, const int& mp, const int& mask) {
        if(mask) {
          auto old=Kokkos::atomic_fetch_add(&counter(0), 1);
          mpAppID_array(old) = mpAppID(mp);
          mpProc_array(old)  = MPs2Proc(mp);
        }
      };
      mps_g->parallel_for(setMPProcAppID, "setProcAppID"); 
      //Copy them to CPU loop over them and find which of them need to be migrated so that the data
      //can be prpepared to be sent from CPU
      auto mpAppID_host = Kokkos::create_mirror_view(mpAppID_array);
      auto mpProc_host  = Kokkos::create_mirror_view(mpProc_array);

      Kokkos::deep_copy(mpAppID_host, mpAppID_array);
      Kokkos::deep_copy(mpProc_host, mpProc_array);
     

      //Perepare Send Buffers 
      int self, nProcs;
      MPI_Comm comm = mps_g ->getMPIComm();
      MPI_Comm_rank(comm, &self);
      MPI_Comm_size(comm, &nProcs);

      std::vector<std::vector<int>> sendList(nProcs);   // which index to migrate
      std::vector<std::vector<int>> sendAppID(nProcs);  // and unique identifier mpAppID
      std::vector<std::vector<std::vector<double>>> sendFields(nProcs);

      for (int i = 0; i < numMPs; ++i) {
        int dest = mpProc_host(i);
        if (dest != self) {
          sendList[dest].push_back(i);
          sendAppID[dest].push_back(mpAppID_host(i));
          sendFields[dest].push_back(CPU_fields[i]);
        }
      }

      // Count number of particles being sent and received from each process
      std::vector<int> sendCounts(nProcs, 0);
      for (int r = 0; r < nProcs; ++r) {
        sendCounts[r] = sendList[r].size();
      }
      std::vector<int> recvCounts(nProcs, 0);
      MPI_Alltoall( sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, comm);
      //Migrate them including deleting particles/reducing aary size and adding particles in the new process

 
    }

};

}


