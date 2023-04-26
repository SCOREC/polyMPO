#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

#include "pmt_utils.hpp"
#include "pmt_mesh.hpp"
#include "pmt_materialPoints.hpp"

namespace polyMpmTest{

#define maxMPsPerElm 8

class MPM{
  private:
    Mesh mesh_;
    IntView elm2MaterialPoints_; 
    IntView materialPoints2Elm_;
    
  public:
    MaterialPoints* MPs;
    MPM(Mesh& mesh, MaterialPoints* inMPs, IntView elm2MPs, IntView MPs2Elm):
        mesh_(mesh), MPs(inMPs), elm2MaterialPoints_(elm2MPs), materialPoints2Elm_(MPs2Elm) {
      MPs->rebuild(materialPoints2Elm_);
    };
    ~MPM() {
      delete MPs;
    }

    Mesh getMesh() { return mesh_; }
    IntView getElm2MPs() { return elm2MaterialPoints_; }
    IntView getMPs2Elm() { return materialPoints2Elm_; }

    void setElm2MPs(IntView elm2MPs) { elm2MaterialPoints_ = elm2MPs; }
    void setMPs2Elm(IntView MPs2Elm) { materialPoints2Elm_ = MPs2Elm; }
};

}//namespace polyMpmTest end

#endif

