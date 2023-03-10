#ifndef POLYMPMTEST_MPM_H
#define POLYMPMTEST_MPM_H

namespace polyMpmTest{

class MPM {
  private:
    Mesh* mesh_;
    MaterialPoints* materialPoints_;
    IntView elm2MaterialPoints_; 
  
  public:
    MPM();
    MPM(Mesh* mesh, MaterialPoints* MP);  

    Mesh getMesh();
    MaterialPoints getMPs();
    IntView getElm2MPs();

    void setElm2MPs(IntView elm2MPs);
};

}//namespace polyMpmTest end

#endif

