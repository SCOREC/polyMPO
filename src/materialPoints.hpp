#ifndef POLYMPMTEST_MATERIALPOINTS_H
#define POLYMPMTEST_MATERIALPOINTS_H

#include "utils.hpp"

namespace polyMpmTest{

class MaterialPoints {
  private:
    int count_;
    Vector2View positions_;
    IntView elmIDs_;
    BoolView isActive_;

  public:
    MaterialPoints(){};
    MaterialPoints( int count,
                    Vector2View positions,
                    IntView elmIDs, 
                    BoolView isActive):
                    count_(count),
                    positions_(positions),
                    elmIDs_(elmIDs),
                    isActive_(isActive){};
 

    int getCount() { return count_; }
    Vector2View getPositions() { return positions_; }
    IntView getElmsIDs() { return elmIDs_; }
    BoolView isActive() { return isActive_; }

//TODO:MUTATOR  
    void T2LTracking(Vector2View dx);
    
};

}
#endif

