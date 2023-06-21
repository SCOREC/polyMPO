#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

#include "pmt_utils.hpp"
#include "pmt_mesh.hpp"
#include "pmt_materialPoints.hpp"

namespace polyMpmTest{

#define maxMPsPerElm 8

class MPMesh{
  private:
    Mesh mesh_;
    
  public:
    MaterialPoints* MPs;
    MPMesh(Mesh& mesh, MaterialPoints* inMPs):
        mesh_(mesh), MPs(inMPs) {
    };
    ~MPMesh() {
      delete MPs;
    }

    Mesh getMesh() { return mesh_; }
    //void setAssembly(DoubleView asmReturn) { mesh_.setAssemblyReturn(asmReturn); }
    
    void CVTTrackingEdgeCenterBased(Vector2View dx, const int printFlag = -1);
    void CVTTrackingElmCenterBased(Vector2View dx, const int printFlag = -1);
    void T2LTracking(Vector2View dx, const int printFlag = -1);
};

}//namespace polyMpmTest end

#endif

