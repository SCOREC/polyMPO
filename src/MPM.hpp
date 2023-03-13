#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

#include "utils.hpp"
#include "mesh.hpp"
#include "materialPoints.hpp"
#include "wachpressBasis.hpp"
namespace polyMpmTest{

class MPM{
  private:
    Mesh* mesh_;
    MaterialPoints* materialPoints_;
    IntView elm2MaterialPoints_; 
    //IntVie
    
  public:
    MPM(){};
    MPM(Mesh* mesh, MaterialPoints* MPs):
        mesh_(mesh),
        materialPoints_(MPs){};

    Mesh getMesh() { return mesh_; }
    MaterialPoints getMPs() { return materialPoints_; }
    IntView getElm2MPs() { return elm2MaterialPoints_; }

    void setElm2MPs(IntView elm2MPs) { elm2MaterialPoints_ = elm2MPs; }
};

}//namespace polyMpmTest end

#endif

