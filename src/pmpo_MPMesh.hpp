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
    void communicateFields(const std::vector<std::vector<double>>& fieldData, const int numEntities, const int numEntries, int mode,
                           std::vector<std::vector<int>>& recvIDVec,  std::vector<std::vector<double>>& recvDataVec);
    void communicate_and_take_halo_contributions(const Kokkos::View<double**>& meshField, int nEntities, int numEntries, int mode, int op);
    
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
    void calculateStrain();
};

}//namespace polyMPO end

#endif

