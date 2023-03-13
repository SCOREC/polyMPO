#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

#include "mesh.hpp"
#include "materialPoints.hpp"

namespace polyMpmTest{

class  {
  private:
    Mesh* mesh_;
    MaterialPoints* materialPoints_;
    IntView elm2MaterialPoints_; 
    IntView 
    
  public:
    MPM();
    MPM(Mesh* mesh, MaterialPoints* MPs);  

    Mesh getMesh();
    MaterialPoints getMPs();
    IntView getElm2MPs();

    void setElm2MPs(IntView elm2MPs);
};

}//namespace polyMpmTest end

#endif

