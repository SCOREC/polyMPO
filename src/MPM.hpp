#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

#include "utils.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "wachpressBasis.hpp"

namespace polyMpmTest{

#define maxMPsPerElm 8

class MPM{
  private:
    Mesh mesh_;
    MaterialPoints materialPoints_;
    IntView elm2MaterialPoints_; 
    IntView materialPoints2Elm_;
    
  public:
    MPM(){};
    MPM(Mesh mesh, MaterialPoints MPs, IntView elm2MPs, IntView MPs2Elm):
        mesh_(mesh),
        materialPoints_(MPs),
        elm2MaterialPoints_(elm2MPs),
        materialPoints2Elm_(MPs2Elm){};

    Mesh getMesh() { return mesh_; }
    MaterialPoints getMPs() { return materialPoints_; }
    IntView getElm2MPs() { return elm2MaterialPoints_; }
    IntView getMPs2Elm() { return materialPoints2Elm_; }

    void setElm2MPs(IntView elm2MPs) { elm2MaterialPoints_ = elm2MPs; }
    void setMPs2Elm(IntView MPs2Elm) { materialPoints2Elm_ = MPs2Elm; }
};

}//namespace polyMpmTest end

#endif

