#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

#include "pmt_utils.hpp"
#include "pmt_mesh.hpp"
#include "pmt_materialPoints.hpp"

namespace polyMPO{

#define maxMPsPerElm 8

class MPMesh{
  public:
    Mesh* p_mesh;
    MaterialPoints* p_MPs;
    
    MPMesh(Mesh* inMesh, MaterialPoints* inMPs):
        p_mesh(inMesh), p_MPs(inMPs) {
    };
    ~MPMesh() {
      delete p_mesh;
      delete p_MPs;
    }

    void CVTTrackingEdgeCenterBased(Vec2dView dx);
    void CVTTrackingElmCenterBased(Vec2dView dx);
    void T2LTracking(Vec2dView dx);
};

}//namespace polyMPO end

#endif

