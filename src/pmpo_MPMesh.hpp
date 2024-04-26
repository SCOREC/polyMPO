#ifndef POLYMPO_MPM_H
#define POLYMPO_MPM_H

#include "pmpo_utils.hpp"
#include "pmpo_mesh.hpp"
#include "pmpo_materialPoints.hpp"

namespace polyMPO{

#define maxMPsPerElm 8

class MPMesh{
  public:
    Mesh* p_mesh;
    MaterialPoints* p_MPs;
    std::map<MaterialPointSlice, std::function<void()>> reconstructSlice = std::map<MaterialPointSlice, std::function<void()>>();
    
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
    void push();

    DoubleView assemblyV0();
    template <MaterialPointSlice index>
    DoubleView wtScaAssembly();
    template <MaterialPointSlice index>
    Vec2dView wtVec2Assembly();
    template <MaterialPointSlice mpfIndex>
    void assembly(bool basisWeightFlag, bool massWeightFlag);

    template<MaterialPointSlice mpSliceIndex>
    void setReconstructSlice();
    void reconstructSlices();

    void printVTP_mesh(int printVTPIndex);
};

}//namespace polyMPO end

#endif

