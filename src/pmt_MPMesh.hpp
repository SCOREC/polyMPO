#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

#include "pmt_utils.hpp"
#include "pmt_mesh.hpp"
#include "pmt_materialPoints.hpp"

namespace polyMpmTest{

#define maxMPsPerElm 8

class MPMesh{
  public:
    Mesh* mesh;
    MaterialPoints* MPs;
    
    MPMesh(Mesh* inMesh, MaterialPoints* inMPs):
        mesh(inMesh), MPs(inMPs) {
    };
    ~MPMesh() {
      delete mesh;
      delete MPs;
    }

    void CVTTrackingEdgeCenterBased(Vector2View dx);
    void CVTTrackingElmCenterBased(Vector2View dx);
    void T2LTracking(Vector2View dx);
};

}//namespace polyMpmTest end

#endif

