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
  private:
   
    bool isPreComputed;
    
  public:
    
    MPMesh() : isPreComputed(false){};
    void computeMatricesAndSolve(); 
    void resetPreComputeFlag();
    Kokkos::View<double*[vec4d_nEntries]> precomputedVtxCoeffs;
    Kokkos::View<double*[vec3d_nEntries][vec4d_nEntries]> precomputedVtxCoeffs_new;

    Mesh* p_mesh;
    MaterialPoints* p_MPs;

    std::map<MeshFieldIndex, std::function<void()>> reconstructSlice = std::map<MeshFieldIndex, std::function<void()>>();
   
    //For MPI Communication
    int numOwnersTot, numHalosTot;
    std::vector<int> numOwnersOnOtherProcs;
    std::vector<int> numHalosOnOtherProcs;
    std::vector<int>haloOwnerProcs;
    
    std::vector<std::vector<int>> haloOwnerLocalIDs;
    
    //Halos send to owners globalIDs which are converted to local IDs of owners
    //Also keep haloLocalID in owner   
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

    void CVTTrackingEdgeCenterBased(Vec2dView dx);
    void CVTTrackingElmCenterBased(const int printVTPIndex = -1);
    void T2LTracking(Vec2dView dx);
    bool push1P();
    void push_ahead();
    void push_swap();
    void push_swap_pos();
    void push();
    void calcBasis(bool use3DArea);

    DoubleView assemblyV0();
    template <MaterialPointSlice index>
    DoubleView wtScaAssembly();
    template <MaterialPointSlice index>
    Vec2dView wtVec2Assembly();
    template <MeshFieldIndex meshFieldIndex>
    void assembly(int order, MeshFieldType type, bool basisWeightFlag, bool massWeightFlag);
    template <MeshFieldIndex meshFieldIndex>
    void assemblyVtx0();
    template <MeshFieldIndex meshFieldIndex>
    void assemblyElm0();
    template <MeshFieldIndex meshFieldIndex>
    void assemblyVtx1();
    template <MeshFieldIndex meshFieldIndex>
    void subAssemblyVtx1(int vtxPerElm, int nCells, int comp, double* array);
    
    void subAssemblyCoeffs(int vtxPerElm, int nCells, double* m11, double* m12, double* m13, double* m14, 
                                                           double* m22, double* m23, double* m24, 
                                                           double* m33, double* m34, 
                                                           double* m44);
    void solveMatrixAndRegularize(int nVertices, double* m11, double* m12, double* m13, double* m14, 
                                                      double* m22, double* m23, double* m24, 
                                                      double* m33, double* m34,
                                                      double* m44);

    template <MeshFieldIndex meshFieldIndex>
    void reconstruct_full();
    void reconstruct_coeff_full();
    void solveMatrix(const Kokkos::View<double**>& vtxMatrices, double& radius, bool scaling);
    void invertMatrix(const Kokkos::View<double**>& vtxMatrices, const double& radius);
 
    template<MeshFieldIndex meshFieldIndex>
    void setReconstructSlice(int order, MeshFieldType type);
    void reconstructSlices();
    
    // Full assembly on GPU
    void assembleField(int vtxPerElm, int nCells, int nVerticesSolve, int nVertices, double* array_sub, double* array_full);
        
    void printVTP_mesh(int printVTPIndex);
};

}//namespace polyMPO end

#endif

